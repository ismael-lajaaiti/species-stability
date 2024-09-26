using GLV
using Random
using DataFrames
using Distributions
using LinearAlgebra
using CairoMakie
include("makie-theme.jl")
set_theme!(theme_minimal())
Random.seed!(123)

S = 30
K_std = 0.3
mu = -1

sigma_weak = 0.1
c_weak = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma_weak / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
)
df_weak, c_weak = get_resistance_resilience(c_weak; alpha = 0.1)

sigma_strong = 1.5
c_strong = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma_strong / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
)
df_strong, c_strong = get_resistance_resilience(c_strong; alpha = 0.1)

# Figure 1 - Plot resistance versus resilience.
w = full_page_width * cm_to_pt
h = w * 0.6 / 2 * width_height_ratio
size = 1.2 .* (w, h)
fig = Figure(; size);
ax = Axis(
    fig[1, 1];
    xlabel = "Resistance",
    ylabel = "Resilience",
    title = "Weak interactions",
    xscale = log10,
    yscale = log10,
)
scatter!(df_weak.resistance, df_weak.resilience; color = df_weak.relative_yield)
ax = Axis(
    fig[1, 2];
    xlabel = "Resistance",
    title = "Strong interactions",
    xscale = log10,
    yscale = log10,
)
scatter!(df_strong.resistance, df_strong.resilience; color = df_strong.relative_yield)
all_relative_yield = vcat(df_strong.relative_yield, df_weak.relative_yield)
limits = extrema(all_relative_yield)
fig[1, 3] = Colorbar(fig[1, 1]; label = "Relative yield", limits)
fig
save_figure("figures/resis-resil-strong-int", fig, size)

# Figure 2 - Plot (A^(-1))_ii elements.
fig = Figure(; size);
bins = 0:0.05:2
ax1 =
    Axis(fig[1, 1]; xlabel = L"-A_{ii}^{-1}", ylabel = "Count", title = "Weak interactions")
v_ii = -Diagonal(inv(c_weak.A)) * ones(richness(c_weak))
hist!(v_ii; bins)
ax2 = Axis(fig[1, 2]; xlabel = L"-A_{ii}^{-1}", title = "Strong interactions")
v_ii = -Diagonal(inv(c_strong.A)) * ones(richness(c_strong))
hist!(v_ii; bins)
save_figure("figures/invA-strong-int", fig, size)
