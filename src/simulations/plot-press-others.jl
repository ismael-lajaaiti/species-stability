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
mu, sigma = 0, 0.3
K_std = 0.3
c = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
)
N_ref = abundance(c)

# Press on focal species, response of others.
kappa = 0.1
Bdiff = zeros(S, S)
tspan = (0, 1_000)
for i in 1:S
    kappa_list = zeros(S)
    kappa_list[i] = kappa
    K_press = (1 .- kappa_list) .* deepcopy(c.K) # Avoid modifying the original K.
    sol = simulate_press(c, K_press, tspan)
    N_press = sol[end]
    delta_N = (N_press .- N_ref) ./ N_ref
    delta_N[i] = 0
    Bdiff[i, :] = delta_N / kappa
end
Bdiff_avg = sum(Bdiff; dims = 1) / (S - 1) |> vec

inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, width / 1.6), fontsize = 8pt);
ax1 =
    Axis(fig[1, 1]; xlabel = "SL", ylabel = "Response to targeted press\non other species")
ry = relative_yield(c)
sorted_indices = sortperm(ry)
ry_min, ry_max = extrema(ry)
ry_val = LinRange(ry_min, ry_max, 100)
scatter!(ry, Bdiff_avg; color = :grey, label = "simulation")
pred = 1 / (S - 1) * (1 ./ ry_val .- 1)
lines!(ry_val, pred; color = :black, label = "prediction")
axislegend()
fig

save("figures/press-others.png", fig)
