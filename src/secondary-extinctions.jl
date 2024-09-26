using AlgebraOfGraphics
using CairoMakie
using DataFrames
using DifferentialEquations
using Distributions
using LinearAlgebra
using LotkaVolterra
using Random
using StatsBase

Random.seed!(1234)

set_theme!(theme_minimal())
function B_mixed(S; mu = -0.1, sigma = 0.12)
    B = rand(Normal(mu, sigma), S, S)
    B[diagind(B)] .= -1
    B
end

S = 50 # Number of species.
mu = -0.15 # Mean interactions.
df = DataFrame(;
    extinctions_true = Float64[], # Number of secondary extinctions.
    extinctions_expected = Float64[],
    eta = Float64[], # Relative yield of the species going extinct.
    resistance = Float64[], # Resistance to press of the species going extinct.
)
create_interaction_matrix(S) = B_mixed(S; mu)
K = rand(Uniform(1, 50), S) # Carrying capacities.
c = create_communities(S, 1; create_interaction_matrix)[S][1]
eta = equilibrium_abundance(c)
B = c.A
A = Diagonal(K) * B * Diagonal(1 ./ K)
V = -inv(A) # Sensitivity matrix.
resistance = eta .* diag(V)
for i in 1:S
    s_expected = V[:, i] .* (K[i] ./ (eta .* K))
    s_expected /= -s_expected[i]
    extinctions_expected = sum(s_expected .<= -0.99)
    eta_0 = deepcopy(eta)
    eta_0[i] = 0
    problem = ODEProblem(lotka_volterra, eta_0, 10_000, c)
    sol = solve(problem)
    eta_end = sol.u[end]
    s_true = (eta_end .- eta) ./ eta
    extinctions_true = sum(s_true .<= -0.99)
    push!(df, (extinctions_true, extinctions_expected, eta[i], resistance[i]))
end

fig = Figure();
ax1 = Axis(
    fig[1, 1];
    xlabel = "Relative yield",
    ylabel = "Number of \n secondary extinctions",
)
scatter!(df.eta, df.extinctions_true)
ax1 = Axis(
    fig[1, 2];
    xlabel = "True number of \n secondary extinctions",
    ylabel = "Expected number of \n secondary extinctions",
)
scatter!(df.extinctions_true, df.extinctions_expected; alpha = 0.5)
x_extrema = extrema(df.extinctions_true) |> collect
lines!(x_extrema, x_extrema; color = :black, label = "y=x")
axislegend(; position = :lt)
label_panels!(fig, 1, 2)
width = full_page_width * cm_to_pt
height = width * 0.7 / width_height_ratio
save_figure(
    "figures/04_secondary-extinctions",
    # "/tmp/plot",
    fig,
    1.2 .* (width, height),
)
fig
