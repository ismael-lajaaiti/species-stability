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

Random.seed!(1234) # Seed for reproducibility.

S = 50 # Number of species.
mu_dict = Dict("Negative" => -0.15, "Negative and positive" => 0.0)
df = DataFrame(;
    resistance = Float64[],
    eta = Float64[],
    reactivity = Float64[],
    K = Float64[],
    received_abs = Float64[],
    received_alg = Float64[],
    interaction = String[],
    s_expected = Float64[],
    s_true = Float64[],
)
for (interaction, mu) in mu_dict
    create_interaction_matrix(S) = B_mixed(S; mu)
    K = rand(Uniform(1, 50), S) # Carrying capacities.
    c = create_communities(S, 1; create_interaction_matrix)[S][1]
    eta = equilibrium_abundance(c)
    B = c.A
    A = Diagonal(K) * B * Diagonal(1 ./ K)
    V = -inv(A) # Sensitivity matrix.
    resistance = eta .* diag(V)
    reactivity = get_reactivity(c)
    s_expected = V[:, 1] .* (K[1] ./ (eta .* K))
    s_expected /= -s_expected[1]
    eta_0 = vcat(0, 1 .+ eta[2:end])
    problem = ODEProblem(lotka_volterra, eta_0, 10000, c)
    sol = solve(problem)
    eta_end = sol.u[end]
    s_true = (eta_end .- eta) ./ eta
    B_nodiag = B - Diagonal(diag(B))
    @info sum(B_nodiag .> 0) / (S * (S - 1))
    B_weighted = B_nodiag # * Diagonal(eta)
    received_abs = sqrt.(vec(sum(B_weighted .^ 2; dims = 2)))
    received_alg = vec(sum(B_weighted; dims = 2))
    append!(
        df,
        DataFrame(;
            resistance,
            eta,
            reactivity,
            K,
            received_abs,
            received_alg,
            interaction,
            s_expected,
            s_true,
        ),
    )
end

# Figure 1. Reactivity against relative yield.
fig = Figure()
interaction = "Negative"
ax1 = Axis(
    fig[1, 1];
    xlabel = "Relative yield",
    ylabel = "Response to pulse",
    title = interaction,
)
df1 = df[df.interaction.==interaction, :]
scatter!(df1.eta, df1.reactivity)
interaction = "Negative and positive"
ax2 = Axis(fig[1, 2]; xlabel = "Relative yield", ylabel = "", title = interaction)
df2 = df[df.interaction.==interaction, :]
scatter!(df2.eta, df2.reactivity)
label_panels!(fig, 1, 2)
width = full_page_width * cm_to_pt
height = width * 0.7 / width_height_ratio
save_figure(
    "figures/01_reactivity-yield",
    # "/tmp/plot",
    fig,
    1.2 .* (width, height),
)

# Figure 2. Resistance against relative yield.
fig = Figure()
interaction = "Negative"
ax1 = Axis(
    fig[1, 1];
    xlabel = "Relative yield",
    ylabel = "Resistance to press",
    title = interaction,
)
df1 = df[df.interaction.==interaction, :]
scatter!(df1.eta, df1.resistance)
interaction = "Negative and positive"
ax2 = Axis(fig[1, 2]; xlabel = "Relative yield", ylabel = "", title = interaction)
df2 = df[df.interaction.==interaction, :]
scatter!(df2.eta, df2.resistance)
label_panels!(fig, 1, 2)
width = full_page_width * cm_to_pt
height = width * 0.7 / width_height_ratio
save_figure(
    "figures/02_resistance-yield",
    # "/tmp/plot",
    fig,
    1.2 .* (width, height),
)

# Figure 3. Expected against true sensitivity to extinctions.
fig = Figure()
interaction = "Negative"
ax1 = Axis(
    fig[1, 1];
    xlabel = "True relative change \n in abundance",
    ylabel = "Expected relative change \n in abundance",
    title = interaction,
)
df1 = df[df.interaction.==interaction, :]
df1_alive = df1[df1.s_true.>-0.99, :]
df1_extinct = df1[df1.s_true.<=-0.99, :]
scatter!(df1_alive.s_true, df1_alive.s_expected; label = "Alive")
scatter!(df1_extinct.s_true, df1_extinct.s_expected; label = "Extinct")
x_extrema = extrema(df1.s_true) |> collect
lines!(x_extrema, x_extrema; color = :black, label = "y=x")
axislegend(; position = :rb)
interaction = "Negative and positive"
ax2 = Axis(
    fig[1, 2];
    xlabel = "True relative change \n in abundance",
    ylabel = "",
    title = interaction,
)
df2 = df[df.interaction.==interaction, :]
df2_alive = df2[df2.s_true.>-0.99, :]
df2_extinct = df2[df2.s_true.<=-0.99, :]
scatter!(df2_alive.s_true, df2_alive.s_expected; label = "Alive")
scatter!(df2_extinct.s_true, df2_extinct.s_expected; label = "Extinct")
x_extrema = extrema(df2.s_true) |> collect
lines!(x_extrema, x_extrema; color = :black, label = "y=x")
label_panels!(fig, 1, 2)
width = full_page_width * cm_to_pt
height = width * 0.7 / width_height_ratio
save_figure(
    "figures/03_secondary-extinctions",
    # "/tmp/plot",
    fig,
    1.2 .* (width, height),
)
fig
