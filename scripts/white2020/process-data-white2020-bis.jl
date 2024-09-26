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

plot_test = 5
df_test = df[df.Plot.==plot_test, :]

"""
Compute the coeffient of variation for a given algae in a given plot.
The coefficient of variation is computed on the de-trended population cover,
that is, on the residuals of the linear regression.
Moreover, the coefficient of variation is computed on the period 5.1 months to 15 months.
Should be applied to unperturbed plots.
"""
function get_cv(df_plot, algae)
    df_tmp = df_plot[5.1 .<= df_plot.Time .<= 15, :]
    data = DataFrame(; time = df_tmp.Time, cover = df_tmp[!, algae])
    ols = lm(@formula(cover ~ time), data)
    std(residuals(ols)) / mean(data.cover)
end
get_cv(df_test, common_algae[2])

function get_cv2(df_plot, df_ctrl, algae)
    ts = unique(df_plot.Time)
    ts[ts.>5]
    lrr = [log_ratio_response(df_plot, df_ctrl, algae, t) for t in ts]
    data = DataFrame(; time = ts, lrr)
    data = data[.!isnan.(data.lrr), :]
    isempty(data) && return NaN
    ols = lm(@formula(lrr ~ time), data)
    cv = std(residuals(ols)) / abs(mean(data.lrr))
end

df_stacked = stack(df, algae_names)
rename!(df_stacked, :variable => :Algae, :value => :Cover)
df_nodist = df_stacked[.!occursin.("Dist", df_stacked.Treatment), :]
df_nodist = combine(groupby(df_nodist, [:Time, :Diversity, :Algae]), :Cover => mean)

function log_ratio_response(df_plot, df_ctrl, algae, t)
    @assert length(unique(df_plot.Diversity)) == 1
    d = unique(df_plot.Diversity)[1] # Plot's diversity.
    F_dist = df_plot[df_plot.Time.==t, algae][1]
    F_dist == 0 && return NaN
    F_ctrl = df_ctrl[
        df_ctrl.Time.==t.&&df_ctrl.Algae.==algae.&&df_ctrl.Diversity.==d,
        :Cover_mean,
    ][1]
    log(F_dist / F_ctrl)
end
log_ratio_response(df_test, df_nodist, common_algae[2], 12)

function get_resilience1(df_plot, df_ctrl, algae)
    ts = unique(df_plot.Time)
    ts[ts.>5]
    lrr = [log_ratio_response(df_plot, df_ctrl, algae, t) for t in ts]
    data = DataFrame(; time = ts, lrr)
    data = data[.!isnan.(data.lrr), :]
    isempty(data) && return NaN
    ols = lm(@formula(lrr ~ time), data)
    resilience = coef(ols)[2]
    lrr5 = log_ratio_response(df_plot, df_ctrl, algae, 5)
    lrr51 = log_ratio_response(df_plot, df_ctrl, algae, 5.1)
    lrr51 > lrr5 && (resilience *= -1)
    resilience
end
resilience = get_resilience1(df_test, df_nodist, common_algae[5])

function get_resilience2(df_plot, df_ctrl, algae)
    ts = unique(df_plot.Time)
    ts[ts.>5]
    lrr = [log_ratio_response(df_plot, df_ctrl, algae, t) for t in ts]
    lrr = lrr[.!isnan.(lrr)]
    isempty(lrr) && return NaN
    @info lrr
    resilience = 1 / mean(lrr)
    resilience
end

function get_resistance(df_plot, df_ctrl, algae)
    lrr_51 = log_ratio_response(df_plot, df_ctrl, algae, 5.1)
    lrr_5 = log_ratio_response(df_plot, df_ctrl, algae, 5)
    lrr_51 - lrr_5
end
resistance = get_resistance(df_test, df_nodist, common_algae[1])

df_cv = DataFrame(; algae = String[], cv = Float64[], diversity = String[])
diversity_list = ["+PGL", "-PGL", "-L", "-G", "-P", "CONTROL"]
for diversity in diversity_list
    plot_ctrl = df[df.Treatment.=="$diversity", "Plot"] |> unique
    for plot in plot_ctrl, algae in common_algae
        cv = get_cv(df[df.Plot.==plot, :], algae)
        push!(df_cv, (algae, cv, diversity))
    end
end
df_cv
df_cv_avg = combine(groupby(df_cv, [:algae, :diversity]), :cv => mean)

df_cv = DataFrame(; algae = String[], cv = Float64[], diversity = String[])
for diversity in diversity_list
    plot_dist = df[df.Treatment.=="$diversity Dist", "Plot"] |> unique
    for plot in plot_dist, algae in common_algae
        cv = get_cv2(df[df.Plot.==plot, :], df_nodist, algae)
        df_cv_avg = combine(groupby(df_cv, [:algae, :diversity]), :cv => mean)
        push!(df_cv, (algae, cv, diversity))
    end
end
df_cv_avg = combine(groupby(df_cv, [:algae, :diversity]), :cv => mean_no_nan => :cv_mean)

mean_no_nan(x) = mean(x[.!isnan.(x)])

df_resistance = DataFrame(; algae = String[], resistance = Float64[], diversity = String[])
for diversity in diversity_list
    plot_dist = df[df.Treatment.=="$diversity Dist", "Plot"] |> unique
    for plot in plot_dist, algae in common_algae
        resistance = get_resistance(df[df.Plot.==plot, :], df_nodist, algae)
        push!(df_resistance, (algae, resistance, diversity))
    end
end
df_resistance_avg = combine(
    groupby(df_resistance, [:algae, :diversity]),
    :resistance => mean_no_nan => :resistance_mean,
)

df_resilience = DataFrame(; algae = String[], resilience = Float64[], diversity = String[])
for diversity in diversity_list
    plot_dist = df[df.Treatment.=="$diversity Dist", "Plot"] |> unique
    for plot in plot_dist, algae in common_algae
        resilience = get_resilience1(df[df.Plot.==plot, :], df_nodist, algae)
        push!(df_resilience, (algae, resilience, diversity))
    end
end
df_resilience
df_resilience_avg = combine(
    groupby(df_resilience, [:algae, :diversity]),
    :resilience => mean_no_nan => :resilience_mean,
)

df_final =
    innerjoin(df_cv_avg, df_resilience_avg, df_resistance_avg; on = [:algae, :diversity])
fig = Figure();
ax = Axis(
    fig[1, 1];
    xlabel = "Coefficient of variation",
    ylabel = "Resilience",
    # yscale = log10,
)
for diversity in unique(df_final.diversity)
    df_tmp = df_final[df_final.diversity.==diversity, :]
    scatter!(df_tmp.resistance_mean, 1 ./ df_tmp.cv_mean; label = diversity)
end
axislegend()
fig
