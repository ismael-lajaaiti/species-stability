using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(12)

# Global parameters.
S = 30
mu_values = LinRange(-0.8, 0.2, 3)
sigma = 0.1
K_std = 0.3
D = LogNormal(log(0.1), 0.5)
kappa_list = rand(D, 100_000)
am = mean(kappa_list)
hm = harmmean(kappa_list)
tspan = (0, 1_000)

# Create figure.
labels = [["A", "B"], ["C", "D"], ["E", "F"]]
inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, 3width / 1.8), fontsize = 8pt);

for (i, mu) in enumerate(mu_values)
    @info mu
    # Create the community.
    c = rand(
        Community,
        S;
        A_ij = Normal(mu / S, sigma / sqrt(S)),
        K_i = Normal(1, K_std),
        interaction = :core,
    )
    N_ref = abundance(c)
    # Press on whole community.
    n_rep = 1_000
    sensitivity_matrix = zeros(n_rep, S)
    for rep in 1:n_rep
        kappa = rand(D, S)
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
        kappa[i] = rand(D)
        K_press = (1 .- kappa) .* deepcopy(c.K) # Avoid modifying the original K.
        sol = simulate_press(c, K_press, tspan)
        N_press = sol[end]
        delta_N = (N_press[i] - N_ref[i]) / N_ref[i]
        sensitivity_matrix[rep, i] = delta_N / -kappa[i]
    end
    sensitivity_sp = vec(mean(sensitivity_matrix; dims = 1))
    # Compute predicted sensitivity to community press.
    ry = relative_yield(c)
    sorted_indices = sortperm(ry)
    ry_min, ry_max = extrema(ry)
    ry_val = LinRange(ry_min, ry_max, 100)
    s_diag_est = 1 ./ ry_val
    pred_stab = s_diag_est .+ am / hm * (1 .- s_diag_est)
    ax1 = Axis(fig[i, 2]; xlabel = "SL")
    scatter!(ry, sensitivity_sp; color = :orangered3, label = "simulation")
    lines!(ry_val, 1 ./ ry_val; color = :black, label = "analytical\nprediction")
    ax1.yreversed = true
    ax2 = Axis(fig[i, 1]; xlabel = "SL", ylabel = "Sensitivity to press (reversed)")
    scatter!(ry, sensitivity_com; color = :goldenrod, label = "simulation")
    lines!(ry_val, pred_stab; color = :black, label = "analytical\nprediction")
    ax2.yreversed = true
    i == 1 && axislegend(; position = :rt)
    l1 = fig[i, 1] = GridLayout()
    l2 = fig[i, 2] = GridLayout()
    Label(fig[i, 3], L"\mu = %$mu"; rotation = 3pi / 2, tellheight = false)
    for (label, layout) in zip(labels[i], [l1, l2])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = label == "C" ? (0, -80, 5, 0) : (0, 5, 5, 0),
            halign = :right,
        )
    end
end

save("figures/simulations/si-sensitivity-mu.png", fig)
# save("figures/simulations/press.svg", fig)
