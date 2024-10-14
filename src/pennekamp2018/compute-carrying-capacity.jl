using Statistics
using CSV
using DataFrames

# Measure species carrying capacities and how they change with temperature.
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
last_day = maximum(df.day)
day_min = last_day - 7 # When do we start to compute sp. carrying capcities.
df = subset(df, :richness => ByRow(==(1)), :day => ByRow(>=(day_min)))

groups = groupby(df, [:predicted_species, :temperature])
df_K = combine(groups, :species_biomass => mean ∘ skipmissing => :K)
rename!(df_K, :predicted_species => :species)

K_relative_diff = []
for row in eachrow(df_K)
    K_ref =
        subset(df_K, :species => ByRow(==(row.species)), :temperature => ByRow(==(15))).K |>
        first
    δK = (row.K - K_ref) / K_ref
    push!(K_relative_diff, δK)
end
df_K.K_relative_diff = K_relative_diff

CSV.write("data/pennekamp2018/carrying-capacity.csv", df_K)
