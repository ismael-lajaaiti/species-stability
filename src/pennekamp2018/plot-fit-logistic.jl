using CSV
using DataFrames
using Statistics
using LsqFit
using CairoMakie
set_theme!(theme_minimal())

# Let's try to fit the curve Loxo at T=15°C.
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
species_list = df.predicted_species |> unique
n_species = length(species_list)
temperature_list = df.temperature |> unique
n_temperature = length(temperature_list)
logistic(t, p) = p[1] ./ (1 .+ (p[1] ./ p[2] .- 1) .* exp.(-p[3] .* t))

day_list = LinRange(1, 37, 100)
fig = Figure(; size = (900, n_species * 150));
ax_list = Any[undef for _ in 1:n_species, _ in 1:n_temperature]
for (i, species) in enumerate(species_list)
    for (j, temperature) in enumerate(temperature_list)
        df_sp =
            subset(
                df,
                :combination => ByRow(==(species)),
                :temperature => ByRow(==(temperature)),
            ) |> dropmissing
        K0 = 0.1 # Order of magnitude of species biomass.
        B0 = mean(subset(df_sp, :day => ByRow(==(1))).species_biomass)
        r0 = 0.1 # Null prior.
        p0 = [K0, B0, r0]
        logistic_fit = curve_fit(logistic, df_sp.day, df_sp.species_biomass, p0)
        fitted_params = coef(logistic_fit)
        ax_list[i, j] =
            ax = Axis(
                fig[i, j];
                xlabel = "Time (days)",
                ylabel = "Species \n biomass (μg/mL)",
                title = i == 1 ? "T=$temperature °C" : "",
                # yscale = log10,
            )
        i == n_species || hidexdecorations!(ax)
        scatter!(ax, df_sp.day, df_sp.species_biomass; label = "data")
        linestyle = any(abs.(fitted_params) .> 100) ? :dash : :solid
        lines!(
            ax,
            day_list,
            logistic(day_list, fitted_params);
            color = :black,
            label = "fit",
            linestyle,
        )
        if j > 1
            hideydecorations!(ax)
            linkyaxes!(ax_list[1], ax)
        end
    end
    Label(fig[i, n_temperature+1], species; rotation = 3pi / 2, tellheight = false)
end
axislegend(; position = :lt)
fig

save("figures/pennekamp2018/fit-logistic.png", fig)

