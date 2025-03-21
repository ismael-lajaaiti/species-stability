using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
using LinearAlgebra
set_theme!(theme_minimal())
Random.seed!(123)

Vii_pred(eta, mu) = (1 + eta * mu) / (1 + mu)
S = 30
sigma = 0.1
mu_val = [-1, -0.5, 0, 0.3]
K_std = 0.3
df = DataFrame(; mu=Float64[], mu0=Float64[], ry=Float64[], Vii=Float64[])
for mu in mu_val
    c = rand(
        Community,
        S;
        A_ij=Normal(mu / S, sigma / sqrt(S)),
        K_i=Normal(1, K_std),
        interaction=:core,
    )
    V = -inv(c.A)
    Vii = V[diagind(V)]
    ry = relative_yield(c)
    a = core_interactions(c)
    mu0 = vec((sum(a, dims=1) .+ 1) ./ S)
    mu = fill(mu, S)
    append!(df, (; mu, mu0, ry, Vii))
end


inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size=(width, width), fontsize=10pt);
ax = Axis(fig[1, 1], xlabel="Relative yield", ylabel=L"V_{00}")
colors = Makie.wong_colors()  # Use a distinguishable color palette
mu_colors = Dict(mu => colors[i] for (i, mu) in enumerate(mu_val))
alpha = 0.7
for mu in mu_val
    color = mu_colors[mu]
    df_mu = subset(df, :mu => ByRow(==(mu)))
    @info mu
    ry_min, ry_max = extrema(df_mu.ry)
    ry_val = LinRange(ry_min, ry_max, 100)
    scatter!(df_mu.ry, df_mu.Vii; label="μS=$mu", color, alpha)
    # scatter!(df_mu.ry, Vii_pred.(df_mu.ry, df_mu.mu0); color, marker=:utriangle)
    lines!(ry_val, Vii_pred.(ry_val, mu / S))
end
fig[0, 1] = axislegend(; orientation=:horizontal)
fig

save("figures/simulations/feedback.pdf", fig)
