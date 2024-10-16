using Statistics
using CSV
using DataFrames
using GLM
using CairoMakie
using StatsBase
set_theme!(theme_minimal())

df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
groups = groupby(df, [:microcosmID, :predicted_species, :temperature])
data = combine(
    groups,
    :species_biomass => log10 ∘ var => :log_variance,
    :species_biomass => log10 ∘ mean => :log_biomass,
)
T_mean = mean(data.temperature)
data.temperature_centered = data.temperature .- T_mean
rename!(data, :predicted_species => :species)

df_K = DataFrame(CSV.File("data/pennekamp2018/carrying-capacity.csv"))
select!(df_K, [:species, :K, :temperature])
transform!(df_K, :K => ByRow(log10) => :log_K)

data = innerjoin(data, df_K; on = [:species, :temperature])

nullmodel = lm(@formula(log_variance ~ log_biomass * temperature), data)
model = lm(@formula(log_variance ~ (log_biomass + log_K) * temperature), data)

aic(model), aic(nullmodel)

x = ftest(nullmodel.model, model.model)

