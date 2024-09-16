import DifferentialEquations: DiscreteCallback, CallbackSet

"""
    time_to_recovery(sol, i, N_eq; alpha = 0.05)

Compute the time to recovery of species `i`.
The time to recovery is defined as the minimum time at which
the species has deviation from its equilibrium abundance
inferior to (1-`alpha`)% of its initial deviation.
Thus the smaller the `alpha`, the longer the time to recovery.
"""
function time_to_recovery(sol, i, N_eq; alpha = 0.05)
    species_trajectory = abs.(Array(sol)[i, :] .- N_eq[i])
    species_trajectory = species_trajectory ./ species_trajectory[1]
    idx = findfirst(species_trajectory .< alpha)
    sol.t[idx]
end

function stochastic_time_to_recovery(
    species_trajectory,
    N_nopress,
    N_press,
    timesteps;
    alpha = 0.05,
)
    species_deviation = abs.(species_trajectory .- N_nopress) / abs(N_nopress - N_press)
    idx = findfirst(species_deviation .< alpha)
    timesteps[idx]
end

"""
    get_resistance_resilience(c::Community; alpha = 0.01, beta = 0.95)

Return a dataframe for with species resistance and resilience for each species of the community.
If community parameters do not allow for full coexistence, the community is assembled,
and the analysis is performed only on the surviving species.

# Keyword argmuents

  - `alpha` is the threshold for `time_to_recovery`
  - `beta` controls the press intensity
"""
function get_resistance_resilience(
    c::Community;
    n_experiments::Int = 1,
    alpha = 0.05,
    beta = 0.05,
    stochastic::Bool = false,
    noise_intensity = 0.05,
    noise!::Function = (du, u, p, t) -> for i in eachindex(du)
        du[i] = noise_intensity * u[i]
    end,
    t_start = 100,
    t_stop = 200,
    t_end = 300,
)
    c = assemble(c)
    S = richness(c)
    @info "Assembled community with $S species."
    N_eq = abundance(c)
    eta = relative_yield(c)
    df = DataFrame(;
        species = Int[],
        resistance = Float64[],
        resilience = Float64[],
        relative_yield = Float64[],
    )
    if !stochastic
        tspan = (0, 10_000)
        for i in 1:S
            K_new = deepcopy(c.K) # Avoid modifying the original K.
            K_new[i] *= 1 - beta
            delta_K = (c.K[i] - K_new[i]) / c.K[i]
            sol_press_on = simulate_press(c, K_new, tspan)
            N_press = sol_press_on[end] # New abundance, under press perturbation.
            resistance = N_eq[i] * delta_K / abs(N_eq[i] - N_press[i])
            N0 = max.(1e-6, N_press) # Allow extinct species to recover.
            sol_press_off = solve(c, N0, tspan; saveat = 1)
            resil = 1 / time_to_recovery(sol_press_off, i, N_eq; alpha)
            push!(df, (i, resistance, resil, eta[i]))
        end
    else
        df.experiment_id = Int[]
        for k in 1:n_experiments
            @info "Running experiment $k"
            for i in 1:S
                sol = simulate_stochastic_press(c, i, noise!; beta, t_start, t_stop, t_end)
                species_trajectory = Array(sol[i, :])
                idx_at_equilibrium = sol.t .< t_start
                idx_under_press = t_start .<= sol.t .< t_stop
                idx_after_press = t_stop .<= sol.t
                N_no_press = mean(species_trajectory[idx_at_equilibrium])
                N_press = mean(species_trajectory[idx_under_press])
                resist = N_no_press * beta / abs(N_no_press - N_press)
                recovery_time = stochastic_time_to_recovery(
                    species_trajectory[idx_after_press],
                    N_no_press,
                    N_press,
                    sol.t[idx_after_press],
                )
                resil = 1 / recovery_time
                push!(df, (i, resist, resil, eta[i], k))
            end
        end
        df = combine(
            groupby(df, :species),
            :resistance => mean,
            :resilience => mean,
            :relative_yield => mean;
            renamecols = false,
        )
    end
    (df, c)
end

"""
    simulate_stochastic_press(
    c::Community,
    i::Int,
    noise!::Function;
    beta = 0.05,
    t_start = 100,
    t_stop = 200,
    t_end = 300,
    )

Simulate a press on species `i` in the community `c`
with a stochastic noise described by the function `noise!`.

# Keywords

  - `beta` press intensity
  - `t_start` time at which the press starts
  - `t_stop` time at which the press stops
  - `t_end` end of the simulation

See also [`get_resistance_resilience`](@ref).
"""
function simulate_stochastic_press(
    c::Community,
    i::Int,
    noise!::Function;
    beta = 0.05,
    t_start = 100,
    t_stop = 200,
    t_end = 300,
)
    K_no_press = deepcopy(c.K)
    K_press = deepcopy(K_no_press)
    K_press[i] = (1 - beta) * K_no_press[i]
    press_on!(integrator) = integrator.p.K = K_press
    press_off!(integrator) = integrator.p.K = K_no_press
    press_start(u, t, integrator) = t in [t_start]
    press_stop(u, t, integrator) = t in [t_stop]
    cb_start = DiscreteCallback(press_start, press_on!)
    cb_stop = DiscreteCallback(press_stop, press_off!)
    callback = CallbackSet(cb_start, cb_stop)
    tstops = [t_start, t_stop]
    solve(c, abundance(c), (0, t_end), noise!; callback, tstops)
end
