using Statistics
using CSV
using DataFrames
using GLM
using CairoMakie
using StatsBase
set_theme!(theme_minimal())

df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
day_start = 7 # Let some time for community to stabilise.
df = subset(df, :day => ByRow(>=(day_start)))
data = combine(
    groupby(df, [:predicted_species, :richness, :temperature]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)

df_K = DataFrame(;
    predicted_species = String[],
    temperature = Float64[],
    K_estimated = Float64[],
)
for gdf in groupby(data, [:predicted_species, :temperature])
    temperature = first(unique(gdf.temperature))
    species = first(unique(gdf.predicted_species))
    sort!(gdf, :richness)
    model = lm(@formula(species_biomass ~ richness), gdf)
    K_estimated = predict(model, gdf)[1]
    push!(df_K, (species, temperature, K_estimated))
end
df_K

CSV.write("data/pennekamp2018/K_linear-model.csv", df_K)

fig = Figure(; size = (800, 450));
T_list = unique(data.temperature)
ax_list = Any[undef for _ in eachindex(T_list)]
for (j, T) in enumerate(T_list)
    ax_list[j] =
        ax = Axis(
            fig[1, j];
            xlabel = "Richness",
            ylabel = "Species biomass",
            title = "T=$(T)°C",
            yscale = log10,
            xticks = 1:6,
        )
    data_T = subset(data, :temperature => ByRow(==(T)))
    for species in unique(data.predicted_species)
        data_sp = subset(data_T, :predicted_species => ByRow(==(species)))
        scatterlines!(data_sp.richness, data_sp.species_biomass; label = species)
    end
    if j > 1
        linkyaxes!(ax_list[1], ax)
        hideydecorations!(ax)
    end
end
fig[0, 1:6] = Legend(fig, ax_list[1]; orientation = :horizontal)
ax_bottom = Axis(
    fig[2, 1:6];
    xlabel = "Temperature (°C)",
    ylabel = "Species \n carrying capacity",
    aspect = AxisAspect(1.8),
    yscale = log10,
)
for gdf in groupby(df_K, :predicted_species)
    scatterlines!(ax_bottom, gdf.temperature, gdf.K_estimated)
end
fig

save("figures/pennekamp2018/biomass-vs-richness.png", fig)
