using Statistics
using CSV
using DataFrames
using CairoMakie
set_theme!(theme_minimal())

df_K = DataFrame(CSV.File("data/pennekamp2018/carrying-capacity.csv"))

fig = Figure();
ax1 = Axis(
    fig[1, 1];
    xlabel = "Temperature (°C)",
    ylabel = "Carrying capacity (log μg/mL)",
    yscale = log10,
)
for species in unique(df_K.species)
    df_sp = subset(df_K, :species => ByRow(==(species)))
    scatterlines!(df_sp.temperature, df_sp.K; label = species)
end
ax2 = Axis(fig[1, 2]; xlabel = "Temperature (°C)", ylabel = "|K - K₁₅| / K₁₅")
for species in unique(df_K.species)
    df_sp = subset(df_K, :species => ByRow(==(species)))
    scatterlines!(df_sp.temperature, abs.(df_sp.K_relative_diff); label = species)
end
fig[0, 1:2] = axislegend("Species"; orientation = :horizontal)
fig
save("figures/pennekamp2018/carrying-capacity.png", fig)
