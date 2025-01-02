using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(1234)

# Create the community.
S = 30
mu, sigma = -0.5, 0.3
K_std = 0.3
c = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
)
N_ref = abundance(c)

function confindence_interval(x; level = 0.9)
    n = length(x)
    std_error = std(x) / sqrt(n)
    t_critical = quantile(TDist(n - 1), 1 - (1 - level) / 2)
    t_critical * std_error
end

# Press on whole community.
D = LogNormal(0, 0.5)
tspan = (0, 1_000)
n_rep = 1_000
sensitivity_matrix = zeros(n_rep, S)
for rep in 1:n_rep
    kappa = rand(D, S)
    K_press = (1 .+ kappa) .* deepcopy(c.K) # Avoid modifying the original K.
    sol = simulate_press(c, K_press, tspan)
    N_press = sol[end]
    delta_N = (N_press .- N_ref) ./ N_ref
    sensitivity_matrix[rep, :] = delta_N ./ kappa
end
sensitivity_com = vec(mean(sensitivity_matrix; dims = 1))
ci_com = [confindence_interval(x) for x in eachcol(sensitivity_matrix)]

# Press on the focal species.
n_rep = 50
sensitivity_matrix = zeros(n_rep, S)
for rep in 1:n_rep, i in 1:S
    kappa = zeros(S)
    kappa[i] = rand(D)
    K_press = (1 .+ kappa) .* deepcopy(c.K) # Avoid modifying the original K.
    sol = simulate_press(c, K_press, tspan)
    N_press = sol[end]
    delta_N = (N_press[i] - N_ref[i]) / N_ref[i]
    sensitivity_matrix[rep, i] = delta_N / kappa[i]
end
sensitivity_sp = vec(mean(sensitivity_matrix; dims = 1))
ci_sp = [confindence_interval(x) for x in eachcol(sensitivity_matrix)]

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

fig = Figure(; size = (600, 300));
ax1 = Axis(
    fig[1, 1];
    xlabel = "Species relative yield",
    ylabel = "Species sensitivity to press",
    title = "Press on focal species",
)
scatter!(ry, sensitivity_sp; color = :coral, label = "simulation")
errorbars!(ry, sensitivity_sp, ci_sp; color = :coral)
# scatter!(ry, s_diag; color = :gray)
lines!(ry_val, 1 ./ ry_val; color = :black, label = "analytical prediction")
axislegend()
ax2 = Axis(fig[1, 2]; xlabel = "Species relative yield", title = "Press on whole community")
scatter!(ry, sensitivity_com; color = :orange2)
errorbars!(ry, sensitivity_com, ci_com; color = :orange2, whiskerwidth = 5)
lines!(ry_val, pred_stab; color = :black)
l1 = fig[1, 1] = GridLayout()
l2 = fig[1, 2] = GridLayout()
for (label, layout) in zip(["A", "B"], [l1, l2])
    Label(
        layout[1, 1, TopLeft()],
        label;
        font = :bold,
        padding = label == "C" ? (0, -80, 5, 0) : (0, 5, 5, 0),
        halign = :right,
    )
end
fig

isdir("figures/simulations") || mkdir("figures/simulations")
save("figures/simulations/press-species-community.png", fig)
save("figures/simulations/press-species-community.svg", fig)
