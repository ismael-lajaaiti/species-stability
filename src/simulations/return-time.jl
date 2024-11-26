using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(1234)

# Create the community.
S = 30
mu, sigma = -0.5, 0.3
K_std = 0.3
c = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
)
N_ref = abundance(c)
ry = relative_yield(c)

function get_tau(sol; alpha = 0.5)
    δN_matrix = (Array(sol) .- N_ref) ./ N_ref
    xi_matrix = δN_matrix ./ N_ref
    xi_matrix ./= xi_matrix[:, 1]
    t_idx = [findfirst(<=(alpha), xi_i) for xi_i in eachrow(xi_matrix)]
    sol.t[t_idx]
end

alpha = 0.1
D = LogNormal(0, 1)
tspan = (0, 100)
n_rep = 500
tau_matrix = zeros(n_rep, S)
for k in 1:n_rep
    xi = rand(D, S)
    x = xi .* ry
    delta_N = x .* c.K
    sol = simulate_pulse(c, delta_N, tspan; saveat = 0.01)
    tau_matrix[k, :] = get_tau(sol; alpha)
end
tau = vec(mean(tau_matrix; dims = 1))

xi_list = rand(D, 100_000)
am = mean(xi_list)
hm = harmmean(xi_list)
a = am / hm
ry_min, ry_max = extrema(ry)
ry_val = LinRange(ry_min, ry_max, 100)
tau_pred = log(1 / alpha) ./ (a .- ry_val .* (a .- 1))

fig = Figure(; size = (500, 400));
ax = Axis(fig[1, 1]; xlabel = "Relative yield", ylabel = "Return time")
scatter!(ry, tau; color = :grey)
lines!(ry_val, tau_pred; color = :black, label = "prediction")
axislegend(; position = :lt)
fig
