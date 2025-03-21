using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
using LinearAlgebra
set_theme!(theme_minimal())
Random.seed!(123)


S = 30
mu, sigma = -1, 0.1
B0_std = 0.003
c = rand(SublinearCommunity, S;
    A_ij=Normal(mu / S, sigma / sqrt(S)),
    B0_i=Normal(0.01, B0_std),
    m_i=Normal(0.4, 0.0),
)
K = carrying_capacity(c)
Beq = abundance(c)

t_start, t_end = 0, 1
saveat = 1
ts = t_start:saveat:t_end
n_ts = length(ts)
x_matrix = zeros(S, n_ts)
xi_0_avg = 0.01
D = LogNormal(log(xi_0_avg), 0.65)
for i in 1:S
    x0 = zeros(S)
    x0[i] = -0.001 * Beq[i]
    sol = simulate_pulse(c, x0, (t_start, t_end); saveat)
    x_matrix[i, :] = (Array(sol)[i, :] .- Beq[i]) ./ abs.(x0[i])
end
resilience = -log.(abs.(x_matrix[:, end])) / t_end ./ (c.r - c.m)

t_start, t_end = 0, 50
saveat = 0.1
ts = t_start:saveat:t_end
n_ts = length(ts)
n_rep = 1_000
x_matrix = zeros(S, n_ts, n_rep)
xi_0_avg = 0.1
D = LogNormal(log(xi_0_avg), 0.65)
for k in 1:n_rep
    x0 = -rand(D, S) .* Beq
    sol = simulate_pulse(c, x0, (t_start, t_end); saveat)
    x_matrix[:, :, k] = (Array(sol) .- Beq) ./ abs.(x0)
end
x_sp_avg = mean(x_matrix, dims=3)
x_sp_avg = dropdims(x_sp_avg, dims=3)
rshort = -log.(abs.(x_sp_avg[:, 2])) / saveat ./ (c.r - c.m)

A = c.A
Aii = relative_selfregulation(c) .* (c.r .- c.m) ./ Beq
a = Diagonal(1 ./ Aii) * (-Diagonal(Aii) + A)
Khat = -a * Beq
ry_hat = Beq ./ Khat
rs = relative_selfregulation(c)
ry = relative_yield(c)

# Press.
c_copy = deepcopy(c)
si_matrix = zeros(S, n_rep)
for k in 1:n_rep
    m_press = (1 .+ rand(D, S)) .* c.m
    c_copy.m = m_press
    Bpress = abundance(c_copy)
    si_matrix[:, k] = (Bpress .- Beq) ./ Beq / (-xi_0_avg)
end
si = mean(si_matrix, dims=2)
si = dropdims(si, dims=2)


inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 17.8cm
fig = Figure(; size=(width, 0.6width), fontsize=8pt);
ga = GridLayout(fig[1, 1])
gb = GridLayout(fig[2, 1])
# gc = GridLayout(fig[3, 1])
ax2 = Axis(ga[1, 1], xlabel="RY", ylabel="Intrinsic recovery rate")
scatter!(ry, resilience)
ablines!(0, 1, color=:black)
ax2 = Axis(ga[1, 2], xlabel="RS", ylabel="")
hideydecorations!(ax2)
scatter!(rs, resilience)
ablines!(0, 1, color=:black)
ax3 = Axis(ga[1, 3], xlabel="RS", ylabel="RY hat")
scatter!(ry, ry_hat)
pred = -(1 .- 1 ./ ry_hat)
colorrange = extrema(pred)
colormap = :algae
ax4 = Axis(gb[1, 1], xlabel="Local net effects", ylabel="Short-term recovery rate")
scatter!(pred, rshort)
ax5 = Axis(gb[1, 2], xlabel="Local net effects", ylabel="Sensitivity to press \n (reversed)")
ax5.yreversed = true
scatter!(pred, si)
# cb = Colorbar(gb[1, 3]; limits=colorrange, colormap, label="Local net effects")
# for (i, row) in enumerate(eachrow(x_sp_avg))
#     lines!(ts, row; color=pred[i], colorrange, colormap)
# end
ax6 = Axis(gb[1, 3], xlabel="RS", ylabel="Local net effects")
scatter!(rs, (1 ./ ry_hat .- 1), label="sublinear")
fig

save("figures/simulations/si-sublinear.svg", fig)

