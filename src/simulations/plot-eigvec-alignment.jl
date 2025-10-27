using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
using LinearAlgebra
set_theme!(theme_minimal())
Random.seed!(12)

# Create the community.
S = 30
mu, sigma = -1, 0.3
K_std = 0.3

inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, width / 1.3), fontsize = 10pt);
ax1 = Axis(fig[1, 1]; xlabel = "SL", ylabel = "Dominant eigenvecor \n alignment, |V1i|")

for i in 1:10
    c = rand(
        Community,
        S;
        A_ij = Normal(mu / S, sigma / sqrt(S)),
        K_i = Normal(1, K_std),
        interaction = :core,
    )
    N_ref = abundance(c)
    ry = relative_yield(c)
    j = Diagonal(ry) * c.A
    eg = eigvecs(j)
    eg_dom = eg[:, end]
    coefs = norm.(eg_dom)
    scatter!(ax1, ry, coefs)
end
fig

save("figures/si-eigvec-alignment.png", fig)
