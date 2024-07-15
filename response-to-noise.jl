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
mu = -0.05 # Mean interactions.
df = DataFrame(;
    amplification = Float64[], # Noise amplification.
    eta = Float64[], # Relative yield of the species going extinct.
)
create_interaction_matrix(S) = B_mixed(S; mu)
c = create_communities(S, 1; create_interaction_matrix)[S][1]
eta = equilibrium_abundance(c)
K = rand(Uniform(1, 50), S) #./ eta # Carrying capacities.
B = c.A
A = Diagonal(K) * B * Diagonal(1 ./ K)
V = -inv(A) # Sensitivity matrix.
resistance = eta .* diag(V)
reactivity = get_reactivity(c)
interactions = sqrt.(sum(inv(B) .^ 2; dims = 2)) |> vec

function glv(u, p, t)
    r, A, K = p
    r .* (1 .+ A * u ./ K) .* u
end

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

fig = Figure();
for (idx, (noise_type, amplitude)) in enumerate(amplitude_dict)
    ax = Axis(
        fig[1, idx];
        xlabel = "Relative yield",
        ylabel = idx == 1 ? "Response to noise (CV)" : "",
        xscale = log10,
        yscale = log10,
        title = noise_type,
    )
    scatter!(eta, amplitude; color = K)
    ax = Axis(
        fig[2, idx];
        xlabel = "Abundance",
        ylabel = idx == 1 ? "Response to noise (CV)" : "",
        yscale = log10,
        xscale = log10,
    )
    scatter!(abundance, amplitude; color = K)
end
Colorbar(fig[1:2, 4]; label = "Carrying capacity", limits = extrema(K))
Label(fig[0, 1:3], "With positive interactions")
width = full_page_width * cm_to_pt
height = width * 0.9 / 2 * width_height_ratio
label_panels!(fig, 2, 3)
save_figure(
    # "figures/05_response-to-noise_mu05",
    "/tmp/plot",
    fig,
    1.2 .* (width, height),
)
fig

