using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
set_theme!(theme_minimal())
Random.seed!(113)

# Create the community.
S = 30
mu, sigma = 0.3, 0.3
K_std = 0.3
n_cons = 5
dependent_sp = 1:n_cons
nondependent_sp = (n_cons+1):S
u = vcat(fill(-0.1, n_cons), fill(0.4, S - n_cons))

function everyone_coexist(S, mu, sigma, u; iter_max = 10_000)
    iter = 1
    N = fill(-1, S)
    while any(N .< 0) && (iter < iter_max)
        global c = rand(
            Community,
            S;
            A_ij = Normal(mu / S, sigma / sqrt(S)),
            K_i = Normal(1, K_std),
            interaction = :core,
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

# Press on whole community.
D = LogNormal(-1, 0.5)
tspan = (0, 1_000)
n_rep = 1_000
sensitivity_matrix = zeros(n_rep, S)
for rep in 1:n_rep
    kappa = 0.1 * rand(D, S)
    u_press = c.u .+ kappa .* c.u
    sol = simulate_press_u(c, u_press, tspan)
    N_press = sol[end]
    delta_N = (N_press .- N_ref) ./ N_ref
    sensitivity_matrix[rep, :] = delta_N ./ kappa
end
sensitivity_com = vec(mean(sensitivity_matrix; dims = 1))

kappa_list = rand(D, 100_000)
am = mean(kappa_list)
hm = harmmean(kappa_list)
amK = mean(c.K)
hmK = harmmean(c.K)
ratio = am / hm
# ratio2 = ratio * (S - 2n_cons) / S

prediction(inverse_eta, ratio) = inverse_eta * (1 - ratio) + ratio

inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, 0.7width), fontsize = 10pt);
ax = Axis(fig[1, 1]; xlabel = "1 / SL", ylabel = "Sensitivity to press (reversed)")
ax.yreversed = true
# hlines!([1]; color = :grey)
# vlines!([1]; color = :grey)
scatter!(
    1 ./ ry[dependent_sp],
    sensitivity_com[dependent_sp];
    label = "obligate species",
    color = :grey,
)
scatter!(
    1 ./ ry[nondependent_sp],
    sensitivity_com[nondependent_sp];
    label = "non-obligate species",
    color = :black,
)
axislegend(; position = :lt)
inv_eta_min, inv_eta_max = extrema(1 ./ ry)
inv_eta_val = LinRange(inv_eta_min, inv_eta_max, 100)
lines!(inv_eta_val, prediction.(inv_eta_val, ratio); color = :black)
fig

save("figures/si-dependent-species.png", fig)
# save("figures/simulations/dependent-species.svg", fig)
