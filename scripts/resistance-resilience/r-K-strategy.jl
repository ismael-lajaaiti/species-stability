using GLV
using DataFrames
using Distributions
using CairoMakie
include("makie-theme.jl")
set_theme!(theme_minimal())
Random.seed!(1234)

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
c.r = 1 ./ c.K .^ 2 # r-K strategy.
N_eq = abundance(c)

function time_to_recovery(sol, i, N_eq; alpha = 0.05)
    species_trajectory = abs.(Array(sol)[i, :] .- N_eq[i])
    species_trajectory = species_trajectory ./ species_trajectory[1]
    idx = findfirst(species_trajectory .< alpha)
    sol.t[idx]
end

eta = relative_yield(c)
df = DataFrame(;
    resistance = Float64[],
    resilience = Float64[],
    relative_yield = Float64[],
    r = Float64[],
)
tspan = (0, 1_000)
for i in 1:S
    K_new = deepcopy(c.K) # Avoid modifying the original K.
    K_new[i] *= 0.9
    delta_K = (c.K[i] - K_new[i]) / c.K[i]
    sol_press_on = simulate_press(c, K_new, tspan)
    N_press = sol_press_on[end] # New abundance, under press perturbation.
    resistance = N_eq[i] * delta_K / abs(N_eq[i] - N_press[i])
    sol_press_off = solve(c, N_press, tspan; saveat = 1)
    resil = 1 / time_to_recovery(sol_press_off, i, N_eq)
    push!(df, (resistance, resil, eta[i], c.r[i]))
end
df

fig = Figure();
ax = Axis(
    fig[1, 1];
    xlabel = "Resistance",
    ylabel = "Resilience",
    xscale = log10,
    yscale = log10,
    title = "Not corrected for r",
)
scatter!(df.resistance, df.resilience; color = df.relative_yield)
ax = Axis(
    fig[1, 2];
    xlabel = "Resistance",
    xscale = log10,
    yscale = log10,
    title = "Corrected",
)
scatter!(df.resistance, df.resilience ./ df.r; color = df.relative_yield)
limits = extrema(df.relative_yield)
fig[1, 3] = Colorbar(fig[1, 1]; label = "Relative yield", limits)
fig

w = full_page_width * cm_to_pt
h = w * 0.7 / 2 * width_height_ratio
save_figure(
    "figures/rK-strategy",
    # "/tmp/plot",
    fig,
    1.0 .* (w, h),
)
