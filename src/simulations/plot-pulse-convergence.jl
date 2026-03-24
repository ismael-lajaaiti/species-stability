using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using DifferentialEquations
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(123)

# Create the community.
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
N_eq = abundance(c)
ry = relative_yield(c)

# Prepare figure.
inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, 1.3width), fontsize = 8pt);

# Simulation parameters.
saveat = 0.5
t_end = 10
t_start = 0
ts = t_start:saveat:t_end
n_ts = length(ts)

n_rep_vals = [10, 100, 1000]
for (i, n_rep) in enumerate(n_rep_vals)
    g = GridLayout(fig[i, 1])
    xi_matrix = zeros(S, n_ts, n_rep)
    xi_0_avg = 0.01
    D = LogNormal(log(xi_0_avg), 0.65)
    for k in 1:n_rep
        xi_0 = rand(D, S)
        delta_N = -xi_0 .* N_eq
        sol = simulate_pulse(c, delta_N, (t_start, t_end); saveat)
        xi_matrix[:, :, k] = (Array(sol) .- N_eq) ./ (N_eq .* xi_0)
    end
    xi_sp_avg = mean(xi_matrix; dims = 3)
    xi_com_avg = vec(mean(xi_sp_avg; dims = 1))
    xi_list = rand(D, 100_000)
    am = mean(xi_list)
    hm = harmmean(xi_list)
    xi_ratio = am / hm
    xi(eta, t) = xi_ratio * (exp(-t) - exp(-eta * t)) + exp(-eta * t)
    R(eta, t) = eta - log(abs(1 - xi_ratio * (1 - exp(-(1 - eta) * t)))) / t
    R2(eta, t) = -log(abs(xi(eta, t))) / t
    # Axis 2. Species short-term recovery rate.
    ry_min, ry_max = extrema(ry)
    ry_vals = LinRange(ry_min, ry_max, 100)
    ax2 = Axis(g[1, 1]; xlabel = "SL", ylabel = "Recovery rate", title = "Short-term")
    duration = 1.0
    idx = findfirst(==(duration), ts)
    r_simu = -(log.(abs.(xi_sp_avg[:, idx]))) ./ duration
    scatter!(ry, r_simu; color = :goldenrod)
    lines!(ry_vals, R.(ry_vals, duration); color = :black)
    # Axis 3. Species long-term recovery rate.
    ax3 = Axis(g[1, 2]; xlabel = "SL", title = "Long-term")
    i != 1 && (ax2.title = "")
    i != 1 && (ax3.title = "")
    duration = 10
    idx = findfirst(==(duration), ts)
    r_simu = -(log.(abs.(xi_sp_avg[:, idx]))) ./ duration
    scatter!(ry, r_simu; color = :orangered3)
    lines!(ry_vals, R.(ry_vals, duration); color = :black)
    Label(fig[i, 2], "n = $n_rep"; rotation = 3pi / 2, tellheight = false)
end

fig

save("figures/si-pulse-convergence.png", fig)
# save("figures/simulations/pulse.svg", fig)
