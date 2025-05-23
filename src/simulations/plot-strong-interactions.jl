using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(14)

# Create the community.
S_pool = 30
mu, sigma = -2, 2
K_std = 0.3

# Setup figure.
inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 17cm
fig = Figure(; size = (width, 2width / 1.8), fontsize = 8pt);
labels = [["A", "B", "C"], ["D", "E", "F"], ["G", "H", "I"]]

for i in 1:3
    S = 0
    while S < 10
        global c = rand(
            Community,
            S_pool;
            A_ij = Normal(mu / S_pool, sigma / sqrt(S_pool)),
            K_i = Normal(1, K_std),
            interaction = :core,
        )
        c = assemble(c)
        N_ref = abundance(c)
        S = richness(c)
    end
    @info S
    @info N_ref
    # Press on whole community.
    D = LogNormal(log(0.1), 0.5)
    tspan = (0, 1_000)
    n_rep = 1_000
    sensitivity_matrix = zeros(n_rep, S)
    for rep in 1:n_rep
        kappa = 0.1rand(D, S)
        K_press = (1 .- kappa) .* deepcopy(c.K) # Avoid modifying the original K.
        sol = simulate_press(c, K_press, tspan)
        N_press = sol[end]
        delta_N = (N_press .- N_ref) ./ N_ref
        sensitivity_matrix[rep, :] = delta_N ./ -kappa
    end
    sensitivity_com = vec(mean(sensitivity_matrix; dims = 1))
    # Press on the focal species.
    n_rep = 50
    sensitivity_matrix = zeros(n_rep, S)
    for rep in 1:n_rep, i in 1:S
        kappa = zeros(S)
        kappa[i] = 0.1rand(D)
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
    # Plot.
    ax1 = Axis(fig[i, 2]; xlabel = "SL")
    scatter!(ry, sensitivity_sp; color = :orangered3, label = "simulation")
    lines!(ry_val, 1 ./ ry_val; color = :black, label = "analytical\nprediction")
    ax1.yreversed = true
    ax2 = Axis(fig[i, 1]; xlabel = "SL", ylabel = "Sensitivity to press (reversed)")
    scatter!(ry, sensitivity_com; color = :goldenrod, label = "simulation")
    lines!(ry_val, pred_stab; color = :black, label = "analytical\nprediction")
    ax2.yreversed = true
    axislegend(; position = :rt)
    l1 = fig[i, 1] = GridLayout()
    l2 = fig[i, 2] = GridLayout()
    l3 = fig[i, 3] = GridLayout()
    ax3 = Axis(fig[i, 3]; title = "Interaction matrix")
    heatmap!(c.A)
    hidedecorations!(ax3)
    colorrange = extrema(c.A)
    Colorbar(fig[i, 4]; colorrange)
    for (label, layout) in zip(labels[i], [l1, l2, l3[1, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = label == "C" ? (0, -80, 5, 0) : (0, 5, 5, 0),
            halign = :right,
        )
    end
end
fig

save("figures/si-strong-interactions.svg", fig)
