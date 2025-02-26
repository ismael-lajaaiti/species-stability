using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using DifferentialEquations
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(123)

# Create communities.
S = 30
mu_range = LinRange(-1, 0, 10)
sigma = 0.3
K_std = 0.3
coms = [
    rand(
        Community,
        S;
        A_ij=Normal(mu / S, sigma / sqrt(S)),
        K_i=Normal(1, K_std),
        interaction=:core,
    )
    for mu in mu_range]

# Pulse.
n_pulse = 100
saveat = 1
t_start, t_end = 0, 10
t_short = 1
t_long = 10
D = LogNormal(log(0.1), 0.6)

df_pulse = DataFrame(; mu=Float64[], r_short=Float64[], r_long=Float64[])
for (i, mu) in enumerate(mu_range)
    c = coms[i]
    Neq = abundance(c)
    for k in 1:n_pulse
        x_i = -rand(D)
        for i in 1:S
            x = zeros(S)
            x[i] = x_i
            sol = simulate_pulse(c, x, (t_start, t_end); saveat)
            x_traj = vec(Array(sol)[i, :]) .- Neq[i]
            r_short = log(abs(x_traj[1] / x_traj[t_short+1])) / t_short
            r_long = log(abs(x_traj[1] / x_traj[t_long+1])) / t_long
            push!(df_pulse, (mu, r_short, r_long))
        end
    end
end

combine(groupby(df_pulse, :mu), :r_short => mean, :r_long => mean)

N_ref = abundance(c)
