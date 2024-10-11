using Distributions
using CairoMakie
using DataFrames
using GLV
set_theme!(theme_minimal())


"""
Complementary effect in the community `c`.
"""
function complementarity_effect(c::Community)
    ry = relative_yield(c)
    @assert all(ry .> 0) # Ensure that the equilibrium is feasible.
    sum(ry) - 1
end
complementarity_effect(c)

"""
Selection effect of the function defined by the vector of `functional_traits` in the community `c`.
Formally, the selection effect writes SE = S cov(ηᵢ/ ∑ⱼηⱼ, tᵢKᵢ).
By default, the function is assumed to be the community biomass
for which all species contribute equally.
We assume that any community function can be decomposed as ϕ = ∑ᵢtᵢNᵢ
where tᵢ is the functional trait of species i.
"""
function selection_effect(c::Community; functional_traits = fill(1, richness(c)))
    S = richness(c)
    @assert length(functional_traits) == S
    ry = relative_yield(c)
    @assert all(ry .> 0) # Ensure that the equilibrium is feasible.
    ry_hat = ry / sum(ry) # Proportion of relative yield.
    tK = functional_traits .* c.K
    S * cov(ry_hat, tK)
end
selection_effect(c)

S = 30
sigma = 0.2
K_std = 0.2
n = 10
n_perturbations = 500
mu_values = LinRange(-1.2, 0, n)
functional_traits = fill(1, S)
df = DataFrame(;
    community_id = Int[],
    ce = Float64[],
    se = Float64[],
    function_stability = Float64[],
)
for (i, mu) in enumerate(mu_values)
    c = rand(
        Community,
        S;
        A_ij = Normal(mu / S, sigma / sqrt(S)),
        K_i = Normal(1, K_std),
        interaction = :core,
    )
    ce = complementarity_effect(c)
    # se = selection_effect(c)
    ry = relative_yield(c)
    ry_hat = ry / sum(ry)
    N = abundance(c)
    tN = functional_traits .* N
    contribution_to_function = tN / sum(tN)
    se = S * cov(1 ./ ry_hat, contribution_to_function)
    for _ in 1:n_perturbations
        function_stability = 1 / return_time(c, functional_traits)
        push!(df, (i, ce, se, function_stability))
    end
end
df = combine(groupby(df, [:community_id, :ce, :se]), :function_stability => mean)

function return_time(c::Community, functional_traits; alpha = 0.1)
    N_eq = abundance(c)
    x = N_eq .* rand(Normal(0, 0.2), S)
    sol = simulate_pulse(c, x, (0, 100); saveat = 0.1)
    δN_array = Array(sol) .- N_eq
    phi_trajectory = abs.(vec(functional_traits' * δN_array))
    phi_trajectory /= phi_trajectory[1]
    below_threshold = phi_trajectory .< alpha
    time_idx = findfirst(below_threshold)
    sol.t[time_idx]
end
return_time(c, fill(1, S))

scatter(df.se, df.ce)
