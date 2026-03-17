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
c1 = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
    r_i = Uniform(0.5, 2),
)
N_eq = abundance(c1)
ry = relative_yield(c1)

# Pulse.
n_rep = 500
saveat = 0.1
t_end = 10
t_start = 0
ts = t_start:saveat:t_end
n_ts = length(ts)
xi_matrix = zeros(S, n_ts, n_rep)
xi_0_avg = 0.01
D = LogNormal(log(xi_0_avg), 0.65)
for k in 1:n_rep
    xi_0 = rand(D, S)
    delta_N = -xi_0 .* N_eq
    sol = simulate_pulse(c1, delta_N, (t_start, t_end); saveat)
    xi_matrix[:, :, k] = (Array(sol) .- N_eq) ./ (N_eq .* xi_0)
end
xi_sp_avg = mean(xi_matrix; dims = 3)
xi_com_avg = vec(mean(xi_sp_avg; dims = 1))
xi_list = rand(D, 100_000)
am = mean(xi_list)
hm = harmmean(xi_list)
xi_ratio = am / hm
duration = 0.1
idx = findfirst(==(duration), ts)
r_simu = -(log.(abs.(xi_sp_avg[:, idx]))) ./ duration
# Press on whole community.
D = LogNormal(log(0.1), 0.5)
tspan = (0, 1_000)
n_rep = 500
sensitivity_matrix = zeros(n_rep, S)
for rep in 1:n_rep
    kappa = rand(D, S)
    K_press = (1 .- kappa) .* deepcopy(c1.K) # Avoid modifying the original K.
    sol = simulate_press(c1, K_press, tspan)
    N_press = sol[end]
    delta_N = (N_press .- N_eq) ./ N_eq
    sensitivity_matrix[rep, :] = delta_N ./ -kappa
end
sensi_com = vec(mean(sensitivity_matrix; dims = 1))

n_com = 10
r_sp_list = []
r_com_list = []
s_sp_list = []
s_com_list = []
for k in 1:n_com
    mu = rand(Uniform(-1, -0.5))
    c = rand(
        Community,
        S;
        A_ij = Normal(mu / S, sigma / sqrt(S)),
        K_i = Normal(1, K_std),
        interaction = :core,
    )
    N_eq = abundance(c)
    ry = relative_yield(c)
    @info extrema(ry)
    # Pulse.
    n_rep = 100
    saveat = 0.1
    t_end = 1
    t_start = 0
    ts = t_start:saveat:t_end
    n_ts = length(ts)
    xi_matrix = zeros(S, n_ts, n_rep) # Species.
    xi_matrix_com = zeros(n_ts, n_rep) # Entire community.
    xi_0_avg = 0.01
    D = LogNormal(log(xi_0_avg), 0.65)
    for k in 1:n_rep
        xi_0 = rand(D, S)
        delta_N = -xi_0 .* N_eq
        sol = simulate_pulse(c, delta_N, (t_start, t_end); saveat)
        xi_matrix[:, :, k] = (Array(sol) .- N_eq) ./ (N_eq .* xi_0)
        xi_matrix_com[:, k] = (sum(Array(sol); dims = 1) .- sum(N_eq)) ./ sum(N_eq)
        xi_matrix_com[:, k] ./= xi_matrix_com[1, k]
    end
    xi_sp_avg = mean(xi_matrix; dims = 3)
    xi_com_avg = mean(xi_matrix_com; dims = 2)
    xi_list = rand(D, 100_000)
    am = mean(xi_list)
    hm = harmmean(xi_list)
    xi_ratio = am / hm
    duration = 0.5
    idx = findfirst(==(duration), ts)
    r_sp = -(log.(abs.(xi_sp_avg[:, idx]))) ./ duration
    r_com = -(log(abs(xi_com_avg[idx]))) / duration
    @info r_com
    append!(r_sp_list, r_sp)
    push!(r_com_list, r_com)
    # Press.
    D = LogNormal(log(0.1), 0.5)
    tspan = (0, 1_000)
    n_rep = 100
    sensitivity_matrix = zeros(n_rep, S)
    s_com = zeros(n_rep)
    for rep in 1:n_rep
        kappa = rand(D, S)
        K_press = (1 .- kappa) .* deepcopy(c.K) # Avoid modifying the original K.
        sol = simulate_press(c, K_press, tspan)
        N_press = sol[end]
        delta_N = (N_press .- N_eq) ./ N_eq
        sensitivity_matrix[rep, :] = delta_N ./ -kappa
        s_com[rep] = (sum(N_press) - sum(N_eq)) / sum(N_eq) / sum(kappa .* c.K)
    end
    s_sp_avg = mean(sensitivity_matrix; dims = 1)
    s_com_avg = mean(s_com)
    push!(s_com_list, s_com_avg)
    append!(s_sp_list, s_sp_avg)
end

inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 11cm
fig = Figure(; size = (width, width), fontsize = 10pt);
ax1 = Axis(
    fig[1, 1];
    ylabel = "Stability to pulse",
    xlabel = "Stability to pulse",
    title = "Absolute recovery",
)
scatter!(1 .- sensi_com, r_simu; color = ry)
ax2 = Axis(fig[1, 2]; title = "Relative recovery", xlabel = "Stability to press")
scatter!(1 .- sensi_com, r_simu ./ c1.r; color = ry)
Colorbar(fig[1, 3]; limits = extrema(ry), label = "SL")
ax3 = Axis(
    fig[2, 1];
    xlabel = "Stability to press",
    ylabel = "Stability to pulse",
    xscale = Makie.pseudolog10,
    title = "Community level",
)
ax4 = Axis(
    fig[2, 2];
    xlabel = "Stability to press",
    xscale = Makie.pseudolog10,
    title = "Species level",
)
for i in 1:n_com
    scatter!(ax3, 1 - s_com_list[i], r_com_list[i])
    scatter!(ax4, 1 .- s_sp_list[1+S*(i-1):S*i], r_sp_list[1+S*(i-1):S*i]; label = "$i")
end
fig[2, 3] = axislegend(; rowgap = -5)
ax1.xticklabelsvisible = false
ax1.yticklabelsvisible = false
ax2.xticklabelsvisible = false
ax2.yticklabelsvisible = false
ax3.xticklabelsvisible = false
ax3.yticklabelsvisible = false
ax4.xticklabelsvisible = false
ax4.yticklabelsvisible = false
fig

save("figures/blur-dimension.svg", fig)
