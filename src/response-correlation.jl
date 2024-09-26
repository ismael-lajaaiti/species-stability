using AlgebraOfGraphics
using CairoMakie
using DataFrames
using DifferentialEquations
using Distributions
using LinearAlgebra
using LotkaVolterra
using Random
using StatsBase
include("makie-theme.jl")
include("scripts/functions.jl")
set_theme!(theme_minimal())
Random.seed!(1234)


S = 50 # Number of species.
mu = -0.05 # Mean interactions.
df = DataFrame(;
    amplification = Float64[], # Noise amplification.
    eta = Float64[], # Relative yield of the species going extinct.
)
create_interaction_matrix(S) = B_mixed(S; mu)
c = create_communities(S, 1; create_interaction_matrix)[S][1]
eta = equilibrium_abundance(c)
K = rand(Uniform(1, 50), S)

# Compute response to pulse and press.
B = c.A
A = Diagonal(K) * B * Diagonal(1 ./ K)
V = -inv(A) # Sensitivity matrix.
resistance = eta .* diag(V)
reactivity = get_reactivity(c)

# Compute response to noise.
abundance = eta .* K
noise_amplitude = 0.01
alpha_dict = Dict("Immigration" => 0, "Demographic" => 0.5, "Environmental" => 1)
amplitude_dict = Dict() # Where to store the results.
for (noise_type, alpha) in alpha_dict
    w = abundance .^ alpha
    w /= sum(w) # Normalize so that sum(w) = 1.
    white_noise(u, p, t) = noise_amplitude * w
    p = (fill(1, S), A, K)
    pb = SDEProblem(glv, white_noise, eta .* K, (0.0, 10_000.0), p)
    solution = solve(pb)
    colnames = vcat([:t], [Symbol("species_$i") for i in 1:S])
    df = DataFrame(solution, colnames)
    select!(df, Not(:t))
    amplitude = [std(col) / (noise_amplitude * mean(col)) for col in eachcol(df)]
    amplitude_dict[noise_type] = amplitude
end

response_to = Dict(
    "Pulse" => reactivity,
    "Press" => 1 ./ resistance,
    "Immigration \n noise" => amplitude_dict["Immigration"],
    "Demographic \n noise" => amplitude_dict["Demographic"],
    "Environmental \n noise" => amplitude_dict["Environmental"],
)
n = length(response_to)

fig = Figure();
ax_grid = Any[undef for _ in 1:n, _ in 1:n]
for (row, (perturbation_x, response_x)) in enumerate(response_to)
    for (col, (perturbation_y, response_y)) in enumerate(response_to)
        if row >= col
            ax_grid[row, col] = Axis(
                fig[row, col];
                xlabel = row == n ? perturbation_y : "",
                ylabel = col == 1 ? perturbation_x : "",
                yscale = log10,
                xscale = log10,
            )
            scatter!(response_y, response_x; color = eta)
            row != n && hidexdecorations!(ax_grid[row, col])
            col != 1 && hideydecorations!(ax_grid[row, col])
        end
    end
end
Label(fig[0, 1:n], rich("With positive interactions"; font = :bold))
Colorbar(fig[1:n, n+1]; label = "Relative yield", limits = extrema(eta))
w = full_page_width * cm_to_pt
h = w * 0.9 / 2 * width_height_ratio
save_figure(
    "figures/06_correlation-responses_mu05",
    # "/tmp/plot",
    fig,
    1.6 .* (w, h),
)
fig

