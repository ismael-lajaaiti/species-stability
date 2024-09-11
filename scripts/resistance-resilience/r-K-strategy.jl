using GLV
using DataFrames
using Distributions
using CairoMakie
include("makie-theme.jl")
include("scripts/resistance-resilience/functions.jl")
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

df, c = get_resistance_resilience(c)
df.r = c.r[df.species]

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
save_figure("figures/rK-strategy", fig, 1.0 .* (w, h))
