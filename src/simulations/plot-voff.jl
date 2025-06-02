using GLV
using CairoMakie
using LinearAlgebra
using DataFrames

S = 30
mu_values = LinRange(-0.1, 0.1, 10)
sigma = 0.1
n_rep = 100

function get_voff(c::Community)
    V = -inv(c.A)
    V_off = V - Diagonal(V)
    sum(V_off) / (S * (S - 1))
end

df = DataFrame(; mu = Float64[], voff = Float64[])
for mu in mu_values, j in 1:n_rep
    c = rand(
        Community,
        S;
        A_ij = Normal(mu / S, sigma / sqrt(S)),
        K_i = Normal(1, K_std),
        interaction = :core,
    )
    voff = get_voff(c)
    push!(df, (mu / S, voff))
end
df_mean = combine(groupby(df, :mu), :voff => mean)

inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, 0.9width), fontsize = 10pt);
ax = Axis(fig[1, 1]; xlabel = "μ/S", ylabel = "Voff")
scatter!(df_mean.mu, df_mean.voff_mean; color = :grey)
ablines!(0, 1; color = :black)
fig

save("figures/voff.png", fig)
