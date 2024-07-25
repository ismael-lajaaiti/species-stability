using CSV
using CairoMakie
using DataFrames
using GLM
using StatsBase

df = DataFrame(CSV.File("data/white2020.csv"))
remove_trailing_whitespace(string) = replace(string, r"\s+$" => "")
new_colnames = remove_trailing_whitespace.(names(df))
rename!(df, new_colnames)

# Find most common algae to have clear trends.
algae_names = names(df)[8:end]
mean_cover = [mean(skipmissing(df[!, algae])) for algae in algae_names]
common_algae = algae_names[sortperm(mean_cover; rev = true)[1:5]] # Keep the top 5.

# Find plots without grazers, with and wihtout perturbations.
plot_nograzer_ctrl = df[df.Treatment.=="+PGL", "Plot"] |> unique
plot_nograzer_dist = df[df.Treatment.=="+PGL Dist", "Plot"] |> unique

remove_zeros(df, colname) = df[df[!, colname].!=0, :]

# Plot time series for the most common algae.
fig = Figure();
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
    "figures/white2020_timeseries-example",
    # "/tmp/plot",
    fig,
    1.2 .* (w, h),
)
fig

function replace_missing!(df)
    for col in algae_names
        df[ismissing.(df[:, col]), col] .= 0.0
    end
    df
end

replace_missing!(df)

fig = Figure();
ax = Axis(
    fig[1, 1];
    xlabel = "Coefficient of variation",
    ylabel = "Response to pulse",
    # yscale = log10,
    # xscale = log10,
)
treatment_list = ["+PGL", "-PGL", "-L", "-G", "-P", "CONTROL"]
for treatment in treatment_list
    # Now we want to relate the coefficient of variation of algea in control plots
    # to their response to the pulse perturbation in the disturbed plots.
    plot_ctrl = df[df.Treatment.=="$treatment", "Plot"] |> unique
    plot_dist = df[df.Treatment.=="$treatment Dist", "Plot"] |> unique
    df_cv = DataFrame(; algae = String[], cv = Float64[], plot = Int[])
    for plot in plot_ctrl
        data = df[df.Plot.==plot, :]
        deleteat!(data, data.Time .< 6)
        for algae in common_algae
            ts = skipmissing(data[!, algae])
            cv = std(ts) / mean(ts)
            push!(df_cv, (algae, cv, plot); promote = true)
        end
    end
    df_cv_avg = combine(groupby(df_cv, :algae), :cv => mean)
    df_pulse = DataFrame(; algae = String[], pulse = Float64[], plot = Int[])
    for plot in plot_dist
        data = df[df.Plot.==plot, :]
        for algae in common_algae
            cover_eq = mean(data[data.Time.<=5, algae])
            data.Y = abs.(data[!, algae] .- cover_eq) / cover_eq
            r_pulse = mean(data[data.Time.>=5.1, :Y])
            push!(df_pulse, (algae, r_pulse, plot))
        end
    end
    df_pulse_avg = combine(groupby(df_pulse, :algae), :pulse => mean)
    df_response = innerjoin(df_cv_avg, df_pulse_avg; on = :algae)
    @info df_response
    scatter!(Float64.(df_response.cv_mean), df_response.pulse_mean; label = treatment)
end
fig[1, 2] = Legend(fig, ax, "Treatments")
w = two_third_page_width * cm_to_pt
h = w / width_height_ratio
save_figure(
    # "figures/white2020_correlation-cv-vs-pulse",
    "/tmp/plot",
    fig,
    1.4 .* (w, h),
)
fig
