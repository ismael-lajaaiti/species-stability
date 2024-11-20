using Statistics
using CSV
using DataFrames
using GLM
using CairoMakie
using StatsBase
set_theme!(theme_minimal())

df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
df_K = DataFrame(CSV.File("data/pennekamp2018/K_linear-model.csv"))
rename!(df_K, :species => :predicted_species)
day_start = 20 # Let some time for community to stabilise.
day_end = 37 # Let some time for community to stabilise.
df = subset(df, :day => ByRow(>=(day_start)))
df = subset(df, :day => ByRow(<=(day_end)))
df.microcosmID |> unique
id = rand(df.microcosmID)
@info id
df = subset(df, :microcosmID => ByRow(==(id)))
temperature = first(unique(df.temperature))
cv(x) = std(skipmissing(x)) / mean(skipmissing(x))^0.72
df_cv = combine(
    groupby(df, :predicted_species),
    :species_biomass => cv,
    [:species_biomass, :day] => reactivity,
    :species_biomass => mean,
)

df_KT = subset(df_K, :temperature => ByRow(==(temperature)))
df_cv = innerjoin(df_cv, df_KT; on = :predicted_species)
df_cv.ry = df_cv.species_biomass_mean ./ df_cv.K_estimated
fig = Figure(; size = (700, 350));
ax1 = Axis(fig[1:2, 1]; xlabel = "Time (days)", ylabel = "Species biomass", yscale = log10)
for species in unique(df.predicted_species)
    df_sp = subset(df, :predicted_species => ByRow(==(species)))
    y = df_sp.species_biomass
    scatterlines!(df_sp.day, y; label = species)
end
ax2 = Axis(fig[1, 2]; xlabel = "Species biomass", ylabel = "Species AR", xscale = log10)
scatter!(df_cv.species_biomass_mean, df_cv.species_biomass_cv)
ax3 = Axis(
    fig[2, 2];
    xlabel = "Species relative yield",
    ylabel = "Species reactivity",
    xscale = log10,
)
scatter!(df_cv.ry, df_cv.species_biomass_cv)
fig[1:2, 3] = Legend(fig, ax1, "Species")
fig

function ar(x)
    x_nomissing = collect(skipmissing(x))
    length(x_nomissing) <= 1 && return NaN
    autocor(x_nomissing, [2])[1]
end

function reactivity(biomass, day; day_end = 4)
    df = DataFrame(; biomass, day)
    df = subset(df, :day => ByRow(<=(day_end)), :biomass => ByRow(>(0)))
    isempty(df) && return missing
    df.log_biomass = log10.(df.biomass)
    model = lm(@formula(log_biomass ~ day), df)
    coef(model)[2]
end

function decorrelation_time(x)
    x_nomissing = collect(skipmissing(x))
    x_nomissing /= mean(x_nomissing)
    length(x_nomissing) <= 1 && return NaN
    ar_vec = autocor(x_nomissing)
    res = findfirst(abs.(ar_vec) .< 0.1)
    return isnothing(res) ? length(ar_vec) : res
end
data = combine(
    groupby(df, [:microcosmID, :predicted_species, :richness]),
    :species_biomass => cv,
    :species_biomass => ar,
    :species_biomass => decorrelation_time,
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
data = subset(data, :species_biomass_cv => ByRow(!isnan))
data = subset(data, :species_biomass_ar => ByRow(!isnan))
data = combine(
    groupby(data, [:predicted_species, :richness]),
    :species_biomass_cv => mean,
    :species_biomass_ar => mean,
    :species_biomass_decorrelation_time => mean => :species_timescale,
    :species_biomass => mean,
)

fig = Figure();
T_list = unique(data.temperature)
ax_list = Any[undef for _ in eachindex(T_list)]
for (j, T) in enumerate(T_list)
    ax = Axis(
        fig[1, j];
        xlabel = "Species biomass",
        ylabel = "Species decorrelation time",
        yscale = log10,
        xscale = log10,
        xticks = 1:6,
    )
    for species in unique(data.predicted_species)
        data_sp = subset(data, :predicted_species => ByRow(==(species)))
        scatterlines!(
            data_sp.species_biomass_mean,
            data_sp.species_timescale;
            label = species,
        )
    end
end
fig[1, 2] = Legend(fig, ax, "Species")
fig

save("figures/pennekamp2018/biomass-vs-richness.png", fig)
