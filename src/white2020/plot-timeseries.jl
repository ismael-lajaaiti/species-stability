using CSV
using DataFrames
using CairoMakie
include("makie-theme.jl")

df = DataFrame(CSV.File("data/white2020_processed.csv"))

df_species = combine(groupby(df, :species), :cover => mean)
df_species = sort(df_species, :cover_mean; rev = true)
n_algae = 5
common_algae = df_species.species[1:n_algae]

plot_disturbance = first(subset(df[df.disturbance, :], :treatment => ByRow(==("pgl"))).plot)
plot_control = first(subset(df[.!df.disturbance, :], :treatment => ByRow(==("pgl"))).plot)

fig = Figure();

ax1 = Axis(fig[1, 1]; xlabel = "Time (Months)", ylabel = "Cover (%)", title = "Disturbance")
df_plot = subset(df, :plot => ByRow(==(plot_disturbance)))
for algae in common_algae
    df_algae = subset(df_plot, :species => ByRow(==(algae)))
    scatter!(df_algae.month, df_algae.cover; label = algae)
    lines!(df_algae.month, df_algae.cover)
end

ax2 = Axis(fig[1, 2]; xlabel = "Time (Months)", ylabel = "Cover (%)", title = "Control")
df_plot = subset(df, :plot => ByRow(==(plot_control)))
for algae in common_algae
    df_algae = subset(df_plot, :species => ByRow(==(algae)))
    scatter!(df_algae.month, df_algae.cover; label = algae)
    lines!(df_algae.month, df_algae.cover)
end

linkyaxes!(ax1, ax2)
hideydecorations!(ax2)




title_vec = ["Disturbed", "Control"]
ax_vec = Vector{Axis}(undef, 2)
for (i, plot_list) in enumerate([plot_nograzer_dist, plot_nograzer_ctrl])
    title = title_vec[i]
    plot_idx = first(plot_list)
    df_plot = df[df.Plot.==plot_idx, :]
    ax_vec[i] = Axis(fig[1, i]; xlabel = "Time (Months)", ylabel = "Cover (%)", title)
    for algae in common_algae
        data = select(df_plot, :Time, algae)
        data.cover = Float64.(data[!, algae])
        scatter!(data.Time, data.cover; label = algae)
        lines!(data.Time, data.cover)
    end
    if i == 1
        vl = vlines!([5.1]; color = :black, linestyle = :dash)
        axislegend(ax_vec[i], [vl], ["Pulse"])
    end
end
linkyaxes!(ax_vec...)
hideydecorations!(ax_vec[2])
fig[1, 3] = Legend(fig, ax_vec[2], "Species")
label_panels!(fig, 1, 2)
w = full_page_width * cm_to_pt
h = w * 0.7 / width_height_ratio
save_figure(
    # "figures/white2020_timeseries-example",
    "/tmp/plot",
    fig,
    1.2 .* (w, h),
)
fig
