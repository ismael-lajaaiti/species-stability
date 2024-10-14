using CSV
using DataFrames
using CairoMakie
set_theme!(theme_minimal())

df = DataFrame(CSV.File("data/pennekamp2018/raw-data.csv"))
df = select(df, Not(:Column1)) # Remove column 1, which contains line numbers.
df.species_biomass = [x == "NA" ? missing : parse(Float64, x) for x in df.species_biomass]

# Sanity checks.
n_experiments = unique(df.microcosmID) |> length # Expected to be 690.
n_species = unique(df.predicted_species) |> length # Expected to be 6.
unique(df.richness)
unique(df.temperature)

CSV.write("data/pennekamp2018/processed-data.csv", df)

