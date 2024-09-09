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


function get_dataframe(S)
    mu_dict = Dict("Negative" => -0.15, "Negative and positive" => 0.0)
    df = DataFrame(;
        eta = Float64[],
        reactivity = Float64[],
        K = Float64[],
        received_abs = Float64[],
        received_alg = Float64[],
        interaction = String[],
    )
    for (interaction, mu) in mu_dict
        create_interaction_matrix(S) = B_mixed(S; mu)
        K = rand(Uniform(1, 50), S) # Carrying capacities.
        c = create_communities(S, 1; create_interaction_matrix)[S][1]
        eta = equilibrium_abundance(c)
        B = c.A
        A = Diagonal(K) * B * Diagonal(1 ./ K)
        reactivity = get_reactivity(c)
        B_nodiag = B - Diagonal(diag(B))
        received_abs = sqrt.(vec(sum(B_nodiag .^ 2; dims = 2)))
        received_alg = vec(sum(B_nodiag; dims = 2))
        append!(
            df,
            DataFrame(; eta, reactivity, K, received_abs, received_alg, interaction),
        )
    end
    df
end

S = 50 # Number of species.
df = get_dataframe(S)

# Figure 1. Net received interactions verses absolute received interactions.
fig = Figure();
interaction = "Negative"
ax1 = Axis(
    fig[1, 1];
    # xlabel = "Relative yield",
    xlabel = L"\sqrt{\sum_{j \neq i} b_{ij}^2}",
    ylabel = L"\sum_{j \neq i} b_{ij}",
    title = interaction,
)
df1 = df[df.interaction.==interaction, :]
scatter!(df1.received_abs, df1.received_alg)
interaction = "Negative and positive"
ax2 = Axis(
    fig[1, 2];
    xlabel = L"\sqrt{\sum_{j \neq i} b_{ij}^2}",
    ylabel = "",
    title = interaction,
)
df2 = df[df.interaction.==interaction, :]
scatter!(df2.received_abs, df2.received_alg)
label_panels!(fig, 1, 2)
width = full_page_width * cm_to_pt
height = width * 0.7 / width_height_ratio
save_figure(
    # "figures/00_interactions-sign",
    "/tmp/plot",
    fig,
    1.2 .* (width, height),
)
fig

# Figure 2. Impact of species richness.
n_rep = 500 # Number of replicates.
S_values = 2:10
df_cor = DataFrame(; S = Int[], correlation = Float64[], interaction = String[])
for S in S_values
    @info "S = $S"
    df_list = []
    for _ in 1:n_rep
        push!(df_list, get_dataframe(S))
    end
    vc = vcat(df_list...)
    df_tmp = combine(
        groupby(vc, :interaction),
        [:received_abs, :received_alg] => cor => :correlation,
    )
    df_tmp.S = fill(S, nrow(df_tmp))
    @info df_tmp
    append!(df_cor, df_tmp)
end

fig = Figure();
ax = Axis(fig[1, 1]; xlabel = "Species richness", ylabel = "Correlation")
for interaction in unique(df_cor.interaction)
    df_tmp = df_cor[df_cor.interaction.==interaction, :]
    scatter!(df_tmp.S, abs.(df_tmp.correlation); label = interaction)
end
fig

fig = Figure();
interaction = "Negative"
ax1 = Axis(
    fig[1, 1];
    # xlabel = "Relative yield",
    xlabel = L"\sqrt{\sum_{j \neq i} b_{ij}^2}",
    ylabel = L"\sum_{j \neq i} b_{ij}",
    title = interaction,
)
df1 = df[df.interaction.==interaction, :]
scatter!(df1.received_abs, df1.received_alg)
interaction = "Negative and positive"
ax2 = Axis(
    fig[1, 2];
    xlabel = L"\sqrt{\sum_{j \neq i} b_{ij}^2}",
    ylabel = "",
    title = interaction,
)
df2 = df[df.interaction.==interaction, :]
scatter!(df2.received_abs, df2.received_alg)
label_panels!(fig, 1, 2)
width = full_page_width * cm_to_pt
height = width * 0.7 / width_height_ratio
save_figure(
    # "figures/00_interactions-sign",
    "/tmp/plot",
    fig,
    1.2 .* (width, height),
)
fig
