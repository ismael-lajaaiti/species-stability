using CSV
using DataFrames
using CairoMakie
set_theme!(theme_minimal())

# Plot a time series of a species in monoculture.
# For example let's consider the microcosmID = 2, which only contains the 'Colp' species.
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
unique(df.predicted_species)

df_colp = subset(df, :combination => ByRow(==("Colp")))
df_colp.microcosmID |> unique
n_temperature = df.temperature |> unique |> length

fig = Figure(; size = (900, 300));
ax_list = Any[undef for _ in 1:n_temperature]
for (j, df_colp_T) in enumerate(groupby(df_colp, :temperature))
    temperature = df_colp_T.temperature |> first |> string
    ax_list[j] =
        ax = Axis(
            fig[1, j];
            xlabel = "Time (days)",
            ylabel = "Species biomass (μg/mL)",
            title = "T = $temperature °C",
            yscale = log10,
        )
    for sdf in groupby(df_colp_T, :microcosmID)
        replicate = sdf.replicate |> first |> string
        scatterlines!(sdf.day, sdf.species_biomass; label = replicate)
    end
    if j > 1
        hideydecorations!(ax)
        linkyaxes!(ax_list[1], ax)
    end
end
fig

save("figures/pennekamp2018/first-ts.png", fig)
