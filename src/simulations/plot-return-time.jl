using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(123)

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
ry = relative_yield(c)

function get_tau(sol, c; alpha = 0.5)
    N_eq = abundance(c)
    δN_matrix = (Array(sol) .- N_eq) ./ N_eq
    xi_matrix = δN_matrix ./ N_eq
    xi_matrix ./= xi_matrix[:, 1]
    t_idx = [findlast(>(alpha), abs.(xi_i)) for xi_i in eachrow(xi_matrix)]
    sol.t[t_idx]
end

function confindence_interval(x; level = 0.9)
    n = length(x)
    std_error = std(x) / sqrt(n)
    t_critical = quantile(TDist(n - 1), 1 - (1 - level) / 2)
    t_critical * std_error
end

function get_return_rate(sol, c, t)
    N_eq = abundance(c)
    δN_matrix = (Array(sol) .- N_eq)
    xi_matrix = δN_matrix ./ N_eq
    xi_matrix ./= xi_matrix[:, 1]
    time_idx = findfirst(==(t), sol.t)
    -log.(abs.(xi_matrix[:, time_idx])) / t
end

# time_values = [0.1, 1, 10]
# title = ["Short-term", "Mid-term", "Long-term"]
time_values = [0.1, 5]
title = ["Short-term", "Long-term"]
D = LogNormal(0, 1)
tspan = (0, 100)
n_rep = 1_000
df_tau = DataFrame(; species = Int64[], ry = Float64[], tau = Float64[], alpha = Float64[])
for k in 1:n_rep
    xi = rand(D, S)
    x = xi .* ry
    delta_N = x .* c.K
    sol = simulate_pulse(c, delta_N, tspan; saveat = 0.01)
    for time in time_values
        return_rate = get_return_rate(sol, c, time)
        append!(df_tau, (; species = 1:S, ry, tau = return_rate, alpha = fill(time, S)))
    end
end
df_tau = combine(
    groupby(df_tau, [:species, :ry, :alpha]),
    :tau => mean,
    :tau => confindence_interval,
)

xi_list = rand(D, 100_000)
am = mean(xi_list)
hm = harmmean(xi_list)
a = am / hm
ry_min, ry_max = extrema(ry)
ry_val = LinRange(ry_min, ry_max, 100)

n_rep = 500
t_start = 0.1
t_end = 10
saveat = 0.1
time_steps = t_start:saveat:t_end
xi_matrix = zeros(S, length(time_steps), n_rep)
N_eq = abundance(c)
for k in 1:n_rep
    xi = rand(D, S)
    x = xi .* ry
    delta_N = x .* c.K
    sol = simulate_pulse(c, delta_N, (t_start, t_end); saveat)
    xi_tmp = (Array(sol) .- N_eq) ./ N_eq
    xi_tmp ./= xi_tmp[:, 1]
    xi_matrix[:, :, k] = xi_tmp
end
xi_mean = mean(xi_matrix; dims = 3)
xi_mean = reshape(xi_mean, size(xi_mean, 1), size(xi_mean, 2))

l = length(time_steps) - 1
xi_diff = zeros(S, l)
for i in 1:S
    for j in 1:l
        xi_diff[i, j] = (xi_mean[i, j] - xi_mean[i, j+1]) / saveat
        abs(xi_diff[i, j]) > 10 && (xi_diff[i, j] = 0)
    end
end
xi_com = vec(mean(xi_mean; dims = 1))

fig = Figure(; size = (600, 550));
g1 = GridLayout(fig[2, 1])
g2 = GridLayout(fig[1, 1])
colormap = :algae
color_list = [:goldenrod3, :indianred]
whiskerwidth = 5
# g3 = GridLayout(fig[3, 1])
ax_list = Any[undef for _ in eachindex(time_values)]
for (i, gdf) in enumerate(groupby(df_tau, :alpha))
    @assert length(unique(gdf.alpha)) == 1
    alpha = first(unique(gdf.alpha))
    color = color_list[i]
    ax_list[i] =
        ax = Axis(
            g1[1, i];
            xlabel = "Species relative yield",
            ylabel = i == 1 ? "Species return rate" : "",
            title = title[i] * " (t = $alpha)",
        )
    scatter!(gdf.ry, gdf.tau_mean; color, label = "simulation")
    errorbars!(gdf.ry, gdf.tau_mean, gdf.tau_confindence_interval; color, whiskerwidth)
end
ax = Axis(
    g2[1, 1];
    xlabel = "Time",
    ylabel = "Distance to equilibrium",
    # title = "Recovery trajectories",
    xscale = log10,
)
for (i, row) in enumerate(eachrow(xi_mean))
    lines!(time_steps, row; color = ry[i], colorrange = extrema(ry), colormap)
end
cb = Colorbar(g2[1, 2]; limits = extrema(ry), label = "Species relative yield", colormap)
# cb.alignmode = Mixed(; right = 0)

# ax = Axis(
#     g2[1, 2];
#     xlabel = "Mean distance to equilibrium",
#     ylabel = "Instantaneous return rate",
#     # title = "Subordination to community",
# )
# poly!(Point2f[(-0.05, -0.1), (-0.05, 0.4), (0.15, 0.4), (0.15, -0.1)]; color = :lightgray)
# for i in 1:S
#     lines!(
#         xi_com[1:end-1],
#         xi_diff[i, 1:end];
#         color = ry[i],
#         colorrange = extrema(ry),
#         colormap,
#     )
# end
# hlines!([0]; color = :black, label = "y = 0")
# # Inset axis.
# a = 10
# inset_ax = Axis(
#     g2[1, 2];
#     width = Relative(0.5),
#     height = Relative(0.5),
#     halign = 0.1,
#     valign = 0.95,
#     backgroundcolor = :lightgray,
# )
# hidedecorations!(inset_ax)
# for i in 1:S
#     lines!(
#         xi_com[a:end-1],
#         xi_diff[i, a:end];
#         color = ry[i],
#         colorrange = extrema(ry),
#         colormap,
#     )
# end
# hlines!([0]; color = :black, label = "y = 0")
# axislegend(; position = :rb)

for (label, layout) in zip(["A", "B", "C"], [g2[1, 1], g1[1, 1], g1[1, 2]])
    Label(
        layout[1, 1, TopLeft()],
        label;
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right,
    )
end
fig

save("figures/simulations/return-time.png", fig)
save("figures/simulations/return-time.svg", fig)
