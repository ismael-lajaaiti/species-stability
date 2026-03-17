using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
set_theme!(theme_minimal())

"""
    create_community(mu, sigma)

Create a community with strong interactions.
Keep only surviving species.
"""
function create_community(mu, sigma; S_pool = 30, S_min = 10)
    S = 0
    while S < S_min
        c = rand(
            Community,
            S_pool;
            A_ij = Normal(mu / S_pool, sigma / sqrt(S_pool)),
            K_i = Normal(1, K_std),
            interaction = :core,
        )
        global c = assemble(c)
        S = richness(c)
    end
    c
end

Random.seed!(12)
# Create communities.
param_list = [(; mu = -1, sigma = 0.3), (; mu = -2, sigma = 2)]
c_list = [create_community(param...) for param in param_list]
# Set up figure.
inch = 96
pt = 4 / 3
cm = inch / 2.54
fig = Figure(; size = (width, width), fontsize = 8pt);
for (k, c) in enumerate(c_list)
    S = richness(c)
    Neq = abundance(c)
    ry = relative_yield(c)
    z_matrix = zeros(S, S) # z[i, j] response of i to the extinction of j.
    tspan = (0, 1_000)
    # Simulate extinctions.
    for i in 1:S
        x = zeros(S)
        x[i] = -Neq[i] # Set species i extinct.
        solution = simulate_pulse(c, x, tspan)
        Neq_new = solution[end]
        z = (Neq_new .- Neq) ./ Neq
        z_matrix[:, i] = z
    end
    z_response = (vec(sum(z_matrix; dims = 2)) .- 1) ./ (S - 1)
    z_impact = (vec(sum(abs.(z_matrix); dims = 1)) .- 1) ./ (S - 1)
    # Plot results.
    ylabel =
        k == 1 ? "Mean response to extinction, z_ext" : "Mean impact of extinction, e_ext"
    ax1 = Axis(fig[k, 1]; xlabel = "SL", ylabel)
    scatter!(ax1, ry, z_response; color = :grey)
    weak_or_strong = k == 1 ? "Weak" : "Strong"
    Label(fig[k, 3], "$weak_or_strong interactions"; tellheight = false, rotation = 3pi / 2)
    ax2 = Axis(fig[k, 2]; xlabel = "SL")
    scatter!(ax2, ry, z_impact; color = :grey)
end
fig

save("figures/si-extinctions.png", fig)
#
