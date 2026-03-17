using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(12)

# Create the community.
S = 30
mu, sigma = -1, 0.3
K_std = 0.3
c = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
)
N_ref = abundance(c)

# Create figure.
inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, width), fontsize = 8pt);
# Press on whole community.
D = LogNormal(log(0.1), 0.5)
tspan = (0, 1_000)
n_rep_vals = [(10, 1), (100, 5), (1_000, 10)]
for (i, n_rep) in enumerate(n_rep_vals)
    n_rep_com, n_rep_sp = n_rep
    # Community-wide press.
    sensitivity_matrix = zeros(n_rep_com, S)
    for rep in 1:n_rep_com
        kappa = rand(D, S)
        K_press = (1 .- kappa) .* deepcopy(c.K) # Avoid modifying the original K.
        sol = simulate_press(c, K_press, tspan)
        N_press = sol[end]
        delta_N = (N_press .- N_ref) ./ N_ref
        sensitivity_matrix[rep, :] = delta_N ./ -kappa
    end
    sensitivity_com = vec(mean(sensitivity_matrix; dims = 1))
    # Targeted press.
    sensitivity_matrix = zeros(n_rep_sp, S)
    for rep in 1:n_rep_sp, i in 1:S
        kappa = zeros(S)
        kappa[i] = rand(D)
        K_press = (1 .- kappa) .* deepcopy(c.K) # Avoid modifying the original K.
        sol = simulate_press(c, K_press, tspan)
        N_press = sol[end]
        delta_N = (N_press[i] - N_ref[i]) / N_ref[i]
        sensitivity_matrix[rep, i] = delta_N / -kappa[i]
    end
    sensitivity_sp = vec(mean(sensitivity_matrix; dims = 1))
    # Compute predicted sensitivity to community press.
    kappa_list = rand(D, 100_000)
    am = mean(kappa_list)
    hm = harmmean(kappa_list)
    ry = relative_yield(c)
    sorted_indices = sortperm(ry)
    ry_min, ry_max = extrema(ry)
    ry_val = LinRange(ry_min, ry_max, 100)
    s_diag_est = 1 ./ ry_val
    pred_stab = s_diag_est .+ am / hm * (1 .- s_diag_est)
    # Plot results.
    ax1 = Axis(fig[i, 2]; xlabel = "SL", title = "Single species (n = $n_rep_sp)")
    scatter!(ry, sensitivity_sp; color = :orangered3, label = "simulation")
    lines!(ry_val, 1 ./ ry_val; color = :black, label = "analytical\nprediction")
    ax1.yreversed = true
    ylabel = i == 2 ? "Sensitivity to press (reversed)" : ""
    ax2 = Axis(fig[i, 1]; xlabel = "SL", ylabel, title = "All species (n = $n_rep_com)")
    scatter!(ry, sensitivity_com; color = :goldenrod, label = "simulation")
    lines!(ry_val, pred_stab; color = :black, label = "analytical\nprediction")
    ax2.yreversed = true
    # axislegend(; position = :rt)
end

fig

save("figures/si-press-convergence.png", fig)
#
