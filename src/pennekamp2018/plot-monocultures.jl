using CSV
using DataFrames
using CairoMakie
set_theme!(theme_minimal())

# Plot a time series of a species in monoculture.
# For example let's consider the microcosmID = 2, which only contains the 'Colp' species.
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
df_mono = subset(df, :richness => ByRow(==(1)))
species_list = df.predicted_species |> unique
n_species = species_list |> length
n_temperature = df.temperature |> unique |> length

fig = Figure(; size = (900, n_species * 150));
ax_list = Any[undef for _ in 1:n_species, _ in 1:n_temperature]
for (i, species) in enumerate(species_list)
    df_species = subset(df_mono, :predicted_species => ByRow(==(species)))
    for (j, df_species) in enumerate(groupby(df_species, :temperature))
        temperature = df_species.temperature |> first |> string
        ax_list[i, j] =
            ax = Axis(
                fig[i, j];
                xlabel = "Time (days)",
                ylabel = "Species \n biomass (μg/mL)",
                title = i == 1 ? "T = $temperature °C" : "",
                yscale = log10,
            )
        for sdf in groupby(df_species, :microcosmID)
            replicate = sdf.replicate |> first |> string
            scatterlines!(sdf.day, sdf.species_biomass; label = replicate)
        end
        i != n_species && hidexdecorations!(ax)
        if j > 1
            hideydecorations!(ax)
            linkyaxes!(ax_list[1], ax)
        end
    end
    Label(fig[i, n_temperature+1], species; rotation = 3pi / 2, tellheight = false)
end
fig

save("figures/pennekamp2018/monocultures.png", fig)
