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

"""
    get_resistance_resilience(c::Community; alpha = 0.01, beta = 0.95)

Return a dataframe for with species resistance and resilience for each species of the community.
If community parameters do not allow for full coexistence, the community is assembled,
and the analysis is performed only on the surviving species.

# Keyword argmuents

  - `alpha` is the threshold for `time_to_recovery`
  - `beta` controls the press intensity
"""
function get_resistance_resilience(c::Community; alpha = 0.05, beta = 0.95)
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
    tspan = (0, 10_000)
    for i in 1:S
        K_new = deepcopy(c.K) # Avoid modifying the original K.
        K_new[i] *= beta
        delta_K = (c.K[i] - K_new[i]) / c.K[i]
        sol_press_on = simulate_press(c, K_new, tspan)
        N_press = sol_press_on[end] # New abundance, under press perturbation.
        resistance = N_eq[i] * delta_K / abs(N_eq[i] - N_press[i])
        N0 = max.(1e-6, N_press) # Allow extinct species to recover.
        sol_press_off = solve(c, N0, tspan; saveat = 1)
        resil = 1 / time_to_recovery(sol_press_off, i, N_eq; alpha)
        push!(df, (i, resistance, resil, eta[i]))
    end
    (df, c)
end
