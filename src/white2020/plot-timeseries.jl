using CSV
using DataFrames
using CairoMakie
using StatsBase
include("makie-theme.jl")

df = DataFrame(CSV.File("data/white2020_processed.csv"))

df_species = combine(groupby(df, :species), :cover => mean)
df_species = sort(df_species, :cover_mean; rev = true)
n_algae = 3
common_algae = df_species.species[1:n_algae]

function process_for_treatment(df, treatment)
    plot_disturbance =
        subset(df[df.disturbance, :], :treatment => ByRow(==(treatment))).plot
    plot_disturbance = unique(plot_disturbance) # Remove duplicates.
    plot_control = subset(df[.!df.disturbance, :], :treatment => ByRow(==(treatment))).plot
    plot_control = unique(plot_control)
    df_treatment = subset(df, :treatment => ByRow(==(treatment)))
    df_treatment = subset(df_treatment, :species => ByRow(in(common_algae)))
    groups = [:month, :disturbance, :species]
    df_treatment = combine(groupby(df_treatment, groups), :cover => mean)
    df_control = df_treatment[.!df_treatment.disturbance, Not(:disturbance)]
    rename!(df_control, :cover_mean => :cover_control)
    df_disturbance = df_treatment[df_treatment.disturbance, Not(:disturbance)]
    rename!(df_disturbance, :cover_mean => :cover_disturbance)
    typical_cover = combine(groupby(df_control, :species), :cover_control => mean)
    on = [:month, :species] # Column on which to do the join.
    df_processed = innerjoin(df_control, df_disturbance; on)
    df_processed = innerjoin(df_processed, typical_cover; on = :species)
    df_processed = transform(
        df_processed,
        [:cover_control, :cover_disturbance, :cover_control_mean] =>
            ByRow((x, y, z) -> abs(x - y) / z) => :cover_difference,
    )
end

size = (600, 4 * 300)
fig = Figure(; size);
treatment_list = unique(df.treatment)
for (i, treatment) in enumerate(treatment_list)
    ax1 = Axis(
        fig[i, 1];
        xlabel = "Time (Months)",
        ylabel = "Cover (%)",
        title = "Disturbance",
    )
    ax2 = Axis(fig[i, 2]; xlabel = "Time (Months)", ylabel = "Cover (%)", title = "Control")
    ax3 = Axis(fig[i, 3]; xlabel = "Time (Months)", ylabel = "Cover difference")
    df_processed = process_for_treatment(df, treatment)
    for algae in unique(df_processed.species)
        df_algae = subset(df_processed, :species => ByRow(==(algae)))
        scatterlines!(ax1, df_algae.month, df_algae.cover_disturbance; label = algae)
        scatterlines!(ax2, df_algae.month, df_algae.cover_control; label = algae)
        scatterlines!(ax3, df_algae.month, df_algae.cover_difference; label = algae)
    end
    linkyaxes!(ax1, ax2)
    hideydecorations!(ax2)
    Label(fig[i, 4], uppercase(treatment); rotation = 3pi / 2, tellheight = false)
end
save("figures/plot-timeseries-example.png", fig)
