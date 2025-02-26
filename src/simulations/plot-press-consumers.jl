using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(123)

# Create the community.
S = 30
mu, sigma = 0.3, 0.3
K_std = 0.3
n_cons = 16
m = -0.1
u = vcat(fill(m, n_cons), fill(1, S - n_cons))

function everyone_coexist(S, mu, sigma, u; iter_max=10_000)
    iter = 1
    N = fill(-1, S)
    while any(N .< 0) && (iter < iter_max)
        global c = rand(
            Community,
            S;
            A_ij=Normal(mu / S, sigma / sqrt(S)),
            K_i=Normal(1, K_std),
            interaction=:core,
            u,
        )
        N = abundance(c)
        iter += 1
    end
    (iter == iter_max) && return nothing
    c
end

c = everyone_coexist(S, mu, sigma, u)
N_ref = abundance(c)
ry = relative_yield(c)

sol = solve(c, rand(S), (0, 1_000))
sol[end]

# Press on whole community.
D = LogNormal(-1, 0.5)
tspan = (0, 1_000)
n_rep = 1_000
sensitivity_matrix = zeros(n_rep, S)
for rep in 1:n_rep
    kappa = rand(D, S)
    # kappa[1:end] .= 0.1
    # kappa[1:n_cons] .= -0.01
    # u_press = c.u .+ kappa .* abs.(c.u) # Avoid modifying the original K.
    u_press = c.u .+ kappa .* c.u # Avoid modifying the original K.
    sol = simulate_press_u(c, u_press, tspan)
    N_press = sol[end]
    delta_N = (N_press .- N_ref) ./ N_ref
    sensitivity_matrix[rep, :] = delta_N ./ kappa
end
sensitivity_com = vec(mean(sensitivity_matrix; dims=1))
# ci_com = [confindence_interval(x) for x in eachcol(sensitivity_matrix)]

V = -inv(c.A)
V * (c.K .* u_press)
i = 1
(N_press[i] - N_ref[i]) / kappa[i]
V[i, i] / ry[i]
delta_u = u_press - u
(sum(V[i, :] .* c.K .* delta_u) - V[i, i] * delta_u[i] * c.K[i]) / N_ref[i]
mean(kappa[2:end]) * (sum(V[i, 2:end] .* c.K[2:end] .* c.u[2:end])) / N_ref[i]
(sum(V[i, 2:end] .* c.K[2:end] .* c.u[2:end] .* kappa[2:end] .* sign.(c.u[2:end]))) / N_ref[i]
V[i, 2:end] .* c.K[2:end] .* c.u[2:end] .* kappa[2:end] .* sign.(c.u[2:end])
(sum(V[i, 2:end] .* c.K[2:end] .* c.u[2:end] .* kappa[2:end] .* 1)) / N_ref[i]
cor(kappa, c.u)
mean(kappa[2:end]) * (sum(V[i, 2:end] .* c.K[2:end] .* c.u[2:end])) / N_ref[i] + V[i, i] * delta_u[i] / ry[i]

mean(kappa[2:end] .* sign.(u[2:end])) * sum(V[i, 2:end] .* c.K[2:end] .* c.u[2:end]) / N_ref[i]
x = kappa[2:end] .* sign.(u[2:end])
y = V[i, 2:end] .* c.K[2:end] .* c.u[2:end]
(mean(kappa[2:end] .* sign.(u[2:end])) * sum(V[i, 2:end] .* c.K[2:end] .* c.u[2:end]) + cov(x, y)) / N_ref[i]


kappa_list = rand(D, 100_000)
am = mean(kappa_list) #* (1 - 2n_cons / S)
hm = harmmean(kappa_list)
amK = mean(c.K)
hmK = harmmean(c.K)
ratio = am / hm
ratio2 = ratio * (S - 2n_cons) / S

pred = []
expected = []
for i in 1:S
    pr = (N_ref[i] - c.u[i] * c.K[i]) / (S * mean(c.u .* c.K))
    exp = mean(V[i, 1:end]) - V[i, i] / S
    push!(pred, pr)
    push!(expected, exp)
end

fig = Figure();
ax = Axis(fig[1, 1])
scatter!(expected, pred, color=c.u)
abline!(0, 1)
fig

function prediction_basal(eta, ratio, u)
    1 + 2 * u * (1 - u / eta) / (mean(c.K) * (S - n_cons * (u + 1)))
end

function prediction(eta, ratio, u)
    u / eta * (sign(u) - ratio) + ratio
end

fig = Figure();
ax = Axis(fig[1, 1], xlabel="Relative yield", ylabel="Sensitivity to press", title="Consumers")
scatter!(relative_yield(c)[1:n_cons], sensitivity_com[1:n_cons])
rymin, rymax = extrema(ry[1:n_cons])
eta_val = LinRange(rymin, rymax, 100)
# lines!(eta_val, prediction.(eta_val, ratio, m), color=:black)
# lines!(eta_val, prediction.(eta_val, ratio2, m), color=:red)
ax = Axis(fig[1, 2], xlabel="Relative yield", title="Basal")
scatter!(relative_yield(c)[n_cons+1:end], sensitivity_com[n_cons+1:end])
rymin, rymax = extrema(ry[n_cons+1:end])
eta_val = LinRange(rymin, rymax, 100)
# lines!(eta_val, prediction.(eta_val, ratio, 1), color=:black)
# lines!(eta_val, prediction_basal.(eta_val, ratio2, -m), color=:red)
fig
