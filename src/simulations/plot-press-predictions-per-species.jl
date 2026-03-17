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
n_rep_vals = [10, 100, 1_000]
for (i, n_rep) in enumerate(n_rep_vals)
    # Community-wide press.
    sensitivity_matrix = zeros(n_rep, S)
    kappa_matrix = zeros(n_rep, S)
    for rep in 1:n_rep
        kappa = rand(D, S)
        kappa_matrix[rep, :] = kappa
        K_press = (1 .- kappa) .* deepcopy(c.K) # Avoid modifying the original K.
        sol = simulate_press(c, K_press, tspan)
        N_press = sol[end]
        delta_N = (N_press .- N_ref) ./ N_ref
        sensitivity_matrix[rep, :] = delta_N ./ -kappa
    end
    sensitivity_com = vec(mean(sensitivity_matrix; dims = 1))
    # Compute predicted sensitivity to community press.
    kappa_list = reduce(vcat, kappa_matrix)
    am = mean(kappa_list)
    hm = harmmean(kappa_list)
    ry = relative_yield(c)
    sorted_indices = sortperm(ry)
    ry_min, ry_max = extrema(ry)
    ry_val = LinRange(ry_min, ry_max, 100)
    s_diag_est = 1 ./ ry
    pred_stab = s_diag_est .+ am / hm * (1 .- s_diag_est)
    # Prediction per species.
    pred_stab_sp = zeros(S)
    for sp in 1:S
        am_sp = (sum(kappa_list) - sum(kappa_matrix[:, sp])) / (S - 1) / n_rep
        hm_sp = harmmean(kappa_matrix[:, sp])
        pred_stab_sp[sp] = 1 / ry[sp] + am_sp / hm_sp * (1 - 1 / ry[sp])
    end
    # Plot results.
    ylabel = i == 2 ? "Predicted sensitivity" : ""
    ax2 = Axis(
        fig[i, 1];
        xlabel = "Observed sensitivity",
        ylabel = "Predicted sensitvity",
        title = "Average across species",
    )
    scatter!(sensitivity_com, pred_stab; color = :grey)
    ablines!(0, 1; color = :black)
    ax1 =
        Axis(fig[i, 2]; xlabel = "Observed sensitivity", ylabel = "", title = "Per species")
    scatter!(sensitivity_com, pred_stab_sp; color = :grey)
    ablines!(0, 1; color = :black)
    Label(fig[i, 3], "n = $n_rep"; tellheight = false)
    if i != 1
        ax1.title = ""
        ax2.title = ""
    end
    if i != 3
        ax1.xlabel = ""
        ax2.xlabel = ""
    end
end

fig

save("figures/si-press-perspecies.png", fig)
#
