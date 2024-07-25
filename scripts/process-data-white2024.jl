using CSV
using CairoMakie
using DataFrames
using GLM
using StatsBase
set_theme!(theme_minimal())

df = DataFrame(CSV.File("data/white2024.csv"))
remove_trailing_whitespace(string) = replace(string, r"\s+$" => "")
new_colnames = remove_trailing_whitespace.(names(df))
rename!(df, new_colnames)

# Find most common algae to have clear trends.
algae_names = names(df)[8:end]
mean_cover = [mean(skipmissing(df[!, algae])) for algae in algae_names]
common_algae = algae_names[sortperm(mean_cover; rev = true)[1:5]] # Keep the top 5.

# Find plots without grazers, with and wihtout perturbations.
plot_nograzer_ctrl = df[df.Treatment.=="A_N0_S0", "Plot"] |> unique
plot_nograzer_n = df[df.Treatment.=="A_N1_S0", "Plot"] |> unique
plot_nograzer_s = df[df.Treatment.=="A_N0_S1", "Plot"] |> unique
plot_nograzer_ns = df[df.Treatment.=="A_N1_S1", "Plot"] |> unique

remove_zeros(df, colname) = df[df[!, colname].!=0, :]

# Plot time series for the most common algae.
fig = Figure();
title_vec = ["Nutrients", "Sediments", "Control"]
ncol = length(title_vec)
ax_vec = Vector{Axis}(undef, ncol)
for (i, plot_list) in enumerate([plot_nograzer_ns, plot_nograzer_s, plot_nograzer_ctrl])
    title = title_vec[i]
    plot_idx = plot_list[1]
    df_plot = df[df.Plot.==plot_idx, :]
    ax_vec[i] = Axis(fig[1, i]; xlabel = "Time (Months)", ylabel = "Cover (%)", title)
    for algae in common_algae
        data = select(df_plot, :Time, algae)
        data.cover = Float64.(data[!, algae])
        scatter!(data.Time, data.cover; label = algae)
        lines!(data.Time, data.cover)
    end
    if i != ncol
        vl_start = vlines!([6]; color = :black, linestyle = :dash)
        vl_end = vlines!([12]; color = :grey, linestyle = :dash)
        if i == 2
            fig[0, 2] = Legend(
                fig,
                [vl_start, vl_end],
                ["Press start", "Press end"];
                tellwidth = false,
                orientation = :horizontal,
            )
        end
    end
end
linkyaxes!(ax_vec...)
for i in 2:ncol
    hideydecorations!(ax_vec[i])
end
fig[1, ncol+1] = Legend(fig, ax_vec[2], "Species")
label_panels!(fig, 1, ncol)
w = full_page_width * cm_to_pt
h = w * 0.7 / width_height_ratio
save_figure(
    "figures/white2024_timeseries-example",
    # "/tmp/plot",
    fig,
    1.2 .* (w, h),
)
fig

treatment = "A_N0_S0"
plot = "6"
df_test = df[df.Plot.==plot, :]

# Compute the coefficient of variation on residuals between months 6 and 12.
function get_cv(df_plot, algae)
    df_press = df_plot[6 .<= df_plot.Time .<= 12, :]
    data = DataFrame(; time = df_press.Time, cover = df_press[!, algae])
    ols = lm(@formula(cover ~ time), data)
    std(residuals(ols)) / mean(data.cover)
end
get_cv(df_test, common_algae[1])

# Add cover algae name as a column and stack the data.
df_stacked = stack(df, algae_names)
rename!(df_stacked, :variable => :Algae, :value => :Cover)
df_ctrl = df_stacked[occursin.("N0_S0", df_stacked.Treatment), :]
df_ctrl = combine(groupby(df_ctrl, [:Time, :Diversity, :Algae]), :Cover => mean)

function log_ratio_response(df_plot, df_ctrl, algae, t)
    @assert length(unique(df_plot.Diversity)) == 1
    d = unique(df_plot.Diversity)[1] # Plot's diversity.
    F_dist = df_plot[df_plot.Time.==t, algae][1]
    F_ctrl = df_ctrl[
        df_ctrl.Time.==t.&&df_ctrl.Algae.==algae.&&df_ctrl.Diversity.==d,
        :Cover_mean,
    ][1]
    log(F_dist / F_ctrl)
end

get_resistance(df_plot, df_ctrl, algae) = log_ratio_response(df_plot, df_ctrl, algae, 12)
get_resistance(df_test, df_ctrl, common_algae[1])

function get_resilience(df_plot, df_ctrl, algae)
    lrr = [log_ratio_response(df_plot, df_ctrl, algae, t) for t in 12:15]
    data = DataFrame(; time = 12:15, lrr)
    ols = lm(@formula(lrr ~ time), data)
    resilience = coef(ols)[2]
end
get_resilience(df_test, df_ctrl, common_algae[1])

data_stability = DataFrame(;
    algae = String[],
    cv = Float64[],
    resistance = Float64[],
    resilience = Float64[],
    treatment = String[],
    plot = String[],
)
for treatment in unique(df.Treatment)
    df_tmp = df[df.Treatment.==treatment, :]
    for plot in unique(df_tmp.Plot)
        df_plot = df_tmp[df_tmp.Plot.==plot, :]
        for algae in common_algae
            cv = get_cv(df_plot, algae)
            resistance = get_resistance(df_plot, df_ctrl, algae)
            resilience = get_resilience(df_plot, df_ctrl, algae)
            if !isnan(resilience) && !isnan(cv)
                push!(data_stability, (algae, cv, resistance, resilience, treatment, plot))
            end
        end
    end
end
df_final = combine(
    groupby(data_stability, [:treatment, :algae]),
    :cv => mean,
    :resistance => mean,
    :resilience => mean,
)
df_final.resilience_mean = df_final.resilience_mean .* sign.(-df_final.resistance_mean)

grazer_list = ["A", "B", "C", "D"]
grazer_dict = Dict(
    "A" => "All removed",
    "B" => "Predator removed",
    "C" => "None removed",
    "D" => "Uncaged control",
)
press_dict = Dict("N1_S0" => "Nutrients", "N0_S1" => "Sediments", "N1_S1" => "Both")
response_short_list = ["cv", "resistance", "resilience"]
response_list = response_short_list .* "_mean"
cartesian = CartesianIndices((1:2, 1:2))
ax_grid = Matrix{Axis}(undef, 2, 2)
fig = Figure();
for (i, (grazer, title)) in enumerate(grazer_dict)
    row, col = Tuple(cartesian[i])
    nested_layout = fig[row, col] = GridLayout()
    for (srow, response1) in enumerate(response_list)
        for (scol, response2) in enumerate(response_list)
            if srow <= scol
                continue
            end
            ax =
                ax_grid[row, col] = Axis(
                    nested_layout[srow-1, scol];
                    xlabel = response_short_list[scol],
                    ylabel = response_short_list[srow],
                )
            for (press, label) in press_dict
                treatment = grazer * "_" * press
                df_tmp = df_final[df_final.treatment.==treatment, :]
                df_tmp = stack(df_tmp, Not([:treatment, :algae]))
                r1 = df_tmp[df_tmp.variable.==response1, :value]
                r2 = df_tmp[df_tmp.variable.==response2, :value]
                scatter!(r2, r1; label)
            end
            srow == 2 && scol == 1 && hidexdecorations!(ax)
            srow == 3 && scol == 2 && hideydecorations!(ax)
        end
    end
    Label(nested_layout[0, 1:2], title; font = :bold)
end
fig[0, 1:2] = Legend(fig, ax_grid[1, 1], "Disturbances"; orientation = :horizontal)
w = full_page_width * cm_to_pt
h = w * 1.3 / width_height_ratio
save_figure(
    "figures/white2024_correlation-responses",
    # "/tmp/plot",
    fig,
    1.3 .* (w, h),
)
fig

