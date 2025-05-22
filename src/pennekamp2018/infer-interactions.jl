using Statistics
using CSV
using DataFrames
using GLM
using CairoMakie
using StatsBase
using Distributions
using MixedModels
set_theme!(theme_minimal())

# Process data.
day_start = 7 # Let some time for community to stabilise.
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
df = subset(df, :day => ByRow(>=(day_start)))
df = subset(
    df,
    :microcosmID => ByRow(!in([49, 229, 275, 327, 353, 359, 696, 261, 312, 406])),
)
# df_K = DataFrame(CSV.File("data/pennekamp2018/K_linear-model.csv"))
df_duo = subset(df, :richness => ByRow(==(2))) # Keep duocultures only.
df_mono = subset(df, :richness => ByRow(==(1))) # Keep monocultures only.
df_duo = combine(
    groupby(df_duo, [:predicted_species, :temperature, :richness, :combination]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
df_mono = combine(
    groupby(df_mono, [:predicted_species, :temperature, :richness, :combination]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
select!(df_mono, [:predicted_species, :temperature], :species_biomass => :B_mono)
df_duo =
    innerjoin(df_duo, df_mono; on = [:predicted_species, :temperature], makeunique = true)
rename!(df_duo, :species_biomass => :B_duo)

species_list = df.predicted_species |> unique |> sort
idx = Dict([sp => i for (i, sp) in enumerate(species_list)])
S = length(species_list)
A = zeros(S, S)
A_normalized = zeros(S, S)
gdf = groupby(df_duo, [:combination, :temperature])
for a in gdf
    a.A = (a.B_duo .- a.B_mono) ./ reverse(a.B_duo)
    a.A_normalized = a.A ./ a.B_mono .* reverse(a.B_mono)
    sp1, sp2 = a.predicted_species
    A[idx[sp1], idx[sp2]] = a.A[1]
    A[idx[sp2], idx[sp1]] = a.A[2]
    A_normalized[idx[sp1], idx[sp2]] = a.A_normalized[1]
    A_normalized[idx[sp2], idx[sp1]] = a.A_normalized[2]
end

# Prepare for saving interaction matrix as a dataframe.
df_A_normalized = DataFrame(A_normalized, :auto)
rename!(df_A_normalized, Symbol.(species_list))
df_A_normalized.row = species_list
select!(df_A_normalized, [:row, Symbol.(species_list)...])
CSV.write("data/pennekamp2018/A_normalized.csv", df_A_normalized)
