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

S = 50 # Number of species.
mu = -0.05 # Mean interactions.
df = DataFrame(;
    s_true = Float64[], # Relative change in species abundance due to extinction.
    s_expected = Float64[], # Expected relative change in species abundance due to extinction.
    eta = Float64[], # Species' relative yield.
)
create_interaction_matrix(S) = B_mixed(S; mu, sigma = 0.10)
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
    eta_0 = deepcopy(eta)
    eta_0[i] = 0
    problem = ODEProblem(lotka_volterra, eta_0, 10_000, c)
    sol = solve(problem)
    eta_end = sol.u[end]
    s_true = (eta_end .- eta) ./ eta
    for j in 1:S
        if i != j
            push!(df, (s_true[j], s_expected[j], eta[j]))
        end
    end
end

df.s_true = abs.(df.s_true)
df.s_expected = abs.(df.s_expected)
data = combine(
    groupby(df, :eta),
    :s_true => mean => :s_true,
    :s_expected => mean => :s_expected,
)
data.abundance = data.eta .* K

fig = Figure();
ax1 = Axis(
    fig[1, 1];
    xlabel = "True sensitivity \n to extinctions",
    ylabel = "Expected sensitivity \n to extinctions",
)
x_extrema = extrema(data.s_true) |> collect
lines!(x_extrema, x_extrema; color = :black, label = "y=x")
scatter!(data.s_true, data.s_expected)
axislegend(; position = :lt)
ax2 = Axis(
    fig[1, 2];
    xlabel = "Relative yield",
    ylabel = "Sensitivity to extinctions",
    yscale = log10,
)
scatter!(data.eta, data.s_true)
ax3 = Axis(
    fig[1, 3];
    xlabel = "Abundance",
    ylabel = "Sensitivity to extinctions",
    yscale = log10,
)
scatter!(data.abundance, data.s_true)
label_panels!(fig, 1, 3)
w = full_page_width * cm_to_pt
h = w * 0.6 / width_height_ratio
save_figure(
    # "figures/04_sensitivity-to-extinctions_mu15-sigma10",
    "/tmp/plot",
    fig,
    1.3 .* (w, h),
)
fig

