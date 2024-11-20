using GLV
using DataFrames
using Distributions
using Random
using CairoMakie
# include("scripts/resistance-resilience/functions.jl")
set_theme!(theme_minimal())
Random.seed!(1234)

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

n_rep = 1_000
a, b = 0.001, 0.002
sd = 0.05
m = 0.5
resistance_mat = zeros(n_rep, S)
for rep in 1:n_rep
    tspan = (0, 1_000)
    beta = rand(Normal(m, sd), S)
    # beta = rand(Uniform(a, b), S)
    K_press = (1 .+ beta) .* deepcopy(c.K) # Avoid modifying the original K.
    sol = simulate_press(c, K_press, tspan)
    N_press = sol[end]
    # delta_N = abs.(N_press .- N_ref) ./ N_ref
    delta_N = (N_press .- N_ref) ./ N_ref
    delta_K = beta
    # resistance_mat[rep, :] = delta_K ./ delta_N
    resistance_mat[rep, :] = delta_N ./ delta_K
end
resistance = vec(mean(resistance_mat; dims = 1))

sensitivity = -inv(c.A) .* c.K' ./ N_ref
s_diag = zeros(S)
s_offdiag = zeros(S)
for i in 1:S
    s_i = sensitivity[i, :]
    s_ii = s_i[i]
    s_diag[i] = s_ii
    s_offdiag[i] = sum(s_i) - s_ii
end

n_rep = 100_000
x = zeros(n_rep, S)
for k in 1:n_rep
    kappa = rand(Normal(m, sd), S)
    for i in 1:S
        res = sum(sensitivity[i, j] * kappa[j] / kappa[i] for j in 1:S if j != i)
        res += sensitivity[i, i]
        x[k, i] = res
    end
end
pred_stab = vec(mean(x; dims = 1))


y = zeros(n_rep)
for k in 1:n_rep
    kappa_1, kappa_2 = rand(Normal(m, sd), 2)
    y[k] = kappa_1 / kappa_2
end
alpha = mean(y)

z = zeros(n_rep, S)
for k in 1:n_rep
    kappa = rand(Normal(m, sd), S)
    kr = kappa[2:end] / kappa[1]
    for i in 1:S
        res = cov(sensitivity[i, collect(j for j in 1:S if j != i)], kr)
        res *= (S - 1)
        z[k, i] = res
    end
end
beta = vec(mean(z; dims = 1))

m_reciprocal = mean(1 ./ rand(Normal(m, sd), 1_000_000))
alpha = m * m_reciprocal
pred_stab = (s_diag .+ alpha * s_offdiag)

ry = relative_yield(c)
sorted_indices = sortperm(ry)
ry_min, ry_max = extrema(ry)
ry_val = LinRange(ry_min, ry_max, 100)
s_diag_est = 1 ./ ry_val
pred_stab = s_diag_est .+ alpha * (1 .- s_diag_est)

fig = Figure(; size = (500, 250));
ax1 = Axis(fig[1, 1]; xlabel = "η", ylabel = "Instability to press")
scatter!(ry, resistance; color = :black)
lines!(ry_val, pred_stab; color = :red)
ax2 = Axis(fig[1, 2]; xlabel = "η", ylabel = "Species sensitivity")
scatter!(ax2, ry, s_diag; label = L"s_i^i")
lines!(ax2, ry_val, 1 ./ ry_val; label = L"\frac{1}{\eta_i}")
scatter!(ax2, ry, s_offdiag; label = L"\sum_{j \neq i} s_j^i")
fig[1, 3] = Legend(fig, ax2)
fig


s = sensitivity[i, collect(j for j in 1:S if j != i)]

cov(rand(Uniform(-0.05, 0.05), S - 1), sensitivity[1, 2:end])


