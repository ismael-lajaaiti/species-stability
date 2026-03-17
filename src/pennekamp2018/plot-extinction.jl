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
day_start = 1 # Let some time for community to stabilise.
day_end = Inf
S_min = 5 # Minimal community richness, for analysis.
df_mono = DataFrame(CSV.File("data/pennekamp2018/df_mono.csv"))
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
df = subset(df, :day => ByRow(>=(day_start)))
df = subset(df, :day => ByRow(<=(day_end)))
df = subset(df, :richness => ByRow(>=(S_min)))
df_avg = combine(
    groupby(df, [:predicted_species, :temperature, :combination, :richness]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
df_avg = innerjoin(df_avg, df_mono; on = [:predicted_species, :temperature])
rename!(df_avg, :B_mono => :K)
rename!(df_avg, :species_biomass => :B)

df_full = subset(df_avg, :richness => ByRow(==(6)))
df_full.sl = df_full.B ./ df_full.K

df_5 = subset(df_avg, :richness => ByRow(==(5)))
species_removed = []
for row in eachrow(df_5)
    combination = split(row.combination, ", ")
    removed = setdiff(species_list, combination) |> first
    push!(species_removed, removed)
end
df_5.species_removed = species_removed


species_list = split(df_full.combination |> first, ", ")
gdf = groupby(df_full, [:predicted_species, :temperature])
for a in gdf
    @info a
    B = a.B |> first
    T = a.temperature |> first
    df_5_T = subset(df_5, :temperature => ByRow(==(T)))
    species = first(a.predicted_species)
    sensitivity = []
    for removed in species_list
        if removed != species
            @info df_removed
            df_removed = subset(
                df_5_T,
                :predicted_species => ByRow(==(species)),
                :species_removed => ByRow(==(removed)),
            )
            B_ext = df_removed.B |> first
            push!(sensitivity, (B - B_ext) / B)
        end
    end
    @info sensitivity
    s_ext = mean(sensitivity)
    a.s_ext = [s_ext]
end
df_full.s_ext = Float64.(df_full.s_ext)

df_full.s_pred = -1 / 5 * (1 ./ df_full.sl .- 1)

scatter(df_full.s_ext, df_full.s_pred; color = df_full.predicted_species)
