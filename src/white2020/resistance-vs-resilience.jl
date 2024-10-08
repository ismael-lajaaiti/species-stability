using CSV
using ColorBrewer
using DataFrames
using CairoMakie
using GLM
using StatsBase
using LinearAlgebra
using CategoricalArrays
using Distributions
include("makie-theme.jl")

df = DataFrame(CSV.File("data/white2020_processed.csv"))

"""
Compute the species resistance under a given treatment.
"""
function get_resistance(df, species, treatment; time = 6.0)
    df_species = subset(
        df,
        :species => ByRow(==(species)),
        :treatment => ByRow(==(treatment)),
        :month => ByRow(==(time)),
    )
    cover_ref = df_species[.!df_species.disturbance, :cover] |> mean
    cover_list = subset(df_species, :disturbance).cover
    cover_list[cover_list.==0] .= NaN
    res = cover_ref ./ abs.(cover_list .- cover_ref)
    [(isnan(x) || isinf(x)) ? missing : x for x in res]
end
get_resistance(df, "cor_ofi", "pg")

function get_resilience(df, species, treatment; timespan = (6, 10))
    df_species = subset(
        df,
        :species => ByRow(==(species)),
        :treatment => ByRow(==(treatment)),
        :month => x -> timespan[1] .<= x .<= timespan[2],
    )
    df_cover_ref = combine(
        groupby(subset(df_species, :disturbance => ByRow(!)), :month),
        :cover => mean => :cover_control,
    )
    groups = groupby(subset(df_species, :disturbance), :plot)
    resilience_list = Array{Float64}(undef, length(groups))
    for (i, df) in enumerate(groups)
        df = innerjoin(df, df_cover_ref; on = :month)
        transform!(
            df,
            [:cover, :cover_control] =>
                ByRow((x, y) -> log10(abs(x - y))) => :log_cover_diff,
        )
        df = subset(df, :log_cover_diff => ByRow(x -> !isinf(x)))
        isempty(df) && return nothing
        resilience = -coef(lm(@formula(log_cover_diff ~ month), df))[2]
        resilience_list[i] = resilience
    end
    resilience_list
end
get_resilience(df, "cor_ofi", "pg")

"""
Compute the species coefficient of variation in control plots.
If `detrended` is true, then the coefficient of variation is the
RMSE of the linear model over the species mean cover.
"""
function cv(df, species, treatment; detrended = true)
    df_species = subset(
        df,
        :species => ByRow(==(species)),
        :treatment => ByRow(==(treatment)),
        :disturbance => ByRow(!),
    )
    groups = groupby(df_species, :plot)
    n_groups = length(groups)
    @info n_groups
    cv_list = Any[undef for _ in 1:n_groups]
    for (i, gdf) in enumerate(groups)
        if !detrended
            cv = std(gdf.cover) / mean(gdf.cover)
        else
            model = lm(@formula(cover ~ month), gdf)
            rmse = sqrt(deviance(model) / nrow(gdf))
            cv = rmse / mean(gdf.cover)
        end
        cv_list[i] = cv
    end
    cv_list
end
cv(df, "cor_ofi", "pg")

df_species = combine(groupby(df, :species), :cover => mean)
df_species = sort(df_species, :cover_mean; rev = true)
n_algae = df.species |> unique |> length # 30
common_algae = df_species.species[1:n_algae]

treatment_list = string.([:pg, :pl, :gl, :pgl, :none])
data =
    DataFrame(; species = String[], treatment = String[], resilience = [], resistance = [])
for treat in treatment_list
    for algae in common_algae
        resistance = get_resistance(df, algae, treat)
        resilience = get_resilience(df, algae, treat)
        isnothing(resilience) && continue
        n = length(resilience)
        m = length(resistance)
        @assert n == m
        treatment = fill(treat, n)
        species = fill(algae, n)
        df_tmp = DataFrame(; species, treatment, resilience, resistance)
        append!(data, df_tmp)
    end
end

data.treatment = categorical(data.treatment)
data.treatment = levels!(data.treatment, ["none", "gl", "pg", "pl", "pgl"])
data = data[.!ismissing.(data.resistance), :]
data.resilience = Float64.(data.resilience)
data.resistance = Float64.(data.resistance)

# Fit the model
model = glm(@formula(resilience ~ resistance * treatment), data, Normal())
predictions = Float64.(predict(model, data; interval = :confidence))
data.resilience_pred = predictions.prediction
data.ci_up = Float64.(predictions.upper)
data.ci_low = Float64.(predictions.lower)

# Fit the overall model (without treatment)
overall_model = glm(@formula(resilience ~ resistance), data, Normal())
predictions = predict(overall_model, data; interval = :confidence)
data.resilience_overall_pred = Float64.(predictions.prediction)
data.overall_ci_up = Float64.(predictions.upper)
data.overall_ci_low = Float64.(predictions.lower)

colors = palette("Set2", 6)

size = (800, 600)
fig = Figure(; size);
ax = Axis(fig[1, 1:5]; xlabel = "Resistance", ylabel = "Resilience")
scatter!(data.resistance, data.resilience; color = (:grey, 0.5))
lines!(
    data.resistance,
    data.resilience_overall_pred;
    color = :black,
    label = "Overall Trend",
    linewidth = 2,
)
data_sorted = sort(data, :resistance)
band!(
    data_sorted.resistance,
    data_sorted.overall_ci_low,
    data_sorted.overall_ci_up;
    color = (:black, 0.3),
    transparency = :true,
)
ax_list = Any[undef for _ in 1:5]
for (i, gdf) in enumerate(groupby(data, :treatment))
    t = gdf.treatment |> unique |> first
    ax_list[i] =
        ax = Axis(fig[2, i]; xlabel = "Resistance", ylabel = "Resilience", title = t)
    @info ax_list
    color = colors[i]
    scatter!(gdf.resistance, gdf.resilience; label = t, alpha = 1, color)
    lines!(gdf.resistance, gdf.resilience_pred; alpha = 0.5, color)
    sort!(gdf, :resistance)
    band!(gdf.resistance, gdf.ci_low, gdf.ci_up; color = (color, 0.5))
    if i > 1
        hideydecorations!(ax)
        linkaxes!(ax, ax_list[1])
    end
end
fig
