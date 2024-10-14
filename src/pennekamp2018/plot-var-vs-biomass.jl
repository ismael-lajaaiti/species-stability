using Statistics
using CSV
using DataFrames
using GLM
using CairoMakie
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

model = lm(@formula(log_variance ~ log_biomass * temperature_centered), data)

fig = Figure(; size = (800, 200));
groups = groupby(data, :temperature)
ax_list = Any[undef for _ in eachindex(groups)]
for (j, gdf) in enumerate(groupby(data, :temperature))
    temperature = gdf.temperature |> first
    ax_list[j] =
        ax = Axis(
            fig[1, j];
            xlabel = "Log species \n biomass",
            ylabel = "Log variance in \n species biomass",
            title = "T = $temperature °C",
        )
    scatter!(gdf.log_biomass, gdf.log_variance)
    x_min, x_max = extrema(skipmissing(gdf.log_biomass))
    x = LinRange(x_min, x_max, 2)
    df_x = DataFrame(; log_biomass = x, temperature_centered = fill(temperature - T_mean, 2))
    y = predict(model, df_x)
    lines!(x, y; color = :black)
    if j > 1
        hideydecorations!(ax)
        linkaxes!(ax_list[1], ax)
    end
end
# fig
save("figures/pennekamp2018/var-vs-biomass.png", fig)
