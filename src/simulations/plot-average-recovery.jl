using GLV
using Statistics
using StatsBase
using DataFrames
using Distributions
using Random
using DifferentialEquations
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(123)

# Set up the figure.
inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, 0.5width), fontsize = 10pt);
n_rep = 500 # Modified, 5_000.
saveat = 0.5
t_end = 10
t_start = 0
ts = t_start:saveat:t_end
n_ts = length(ts)

"""
Our expectation of the trajectory of average species deviation.
"""
expected_recovery(t) = -exp(-t)

S = 30
sigma = 0.3
mu_vals = [-1, -0.5, 0]
K_std = 0.3
n_rep = 500
for (i, mu) in enumerate(mu_vals)
    # Create the community.
    c = rand(
        Community,
        S;
        A_ij = Normal(mu / S, sigma / sqrt(S)),
        K_i = Normal(1, K_std),
        interaction = :core,
    )
    N_eq = abundance(c)
    # Simulate the dynamics.
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
    # Plot the recovery.
    ax = Axis(fig[1, i]; xlabel = "Time", ylabel = "Community deviation, <z>")
    ax.title = "mu = $mu"
    scatter!(ax, ts, xi_com_avg; color = :grey)
    lines!(ax, ts, expected_recovery.(ts); color = :black)
    i == 1 || hideydecorations!(ax)
end

fig

save("figures/si-average-recovery.png", fig)
