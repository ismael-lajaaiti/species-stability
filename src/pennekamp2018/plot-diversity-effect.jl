using Statistics
using CSV
using DataFrames
using GLM
using CairoMakie
using StatsBase
using Distributions
using MixedModels
set_theme!(theme_minimal())

# Process data.
day_start = 7 # Let some time for community to stabilise.
S_min = 1 # Minimal community richness, for analysis.
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
df = subset(df, :day => ByRow(>=(day_start)))
df = subset(df, :richness => ByRow(>=(S_min)))
df = subset(
    df,
    :microcosmID => ByRow(!in([49, 229, 275, 327, 353, 359, 696, 261, 312, 406])),
)
df_K = DataFrame(CSV.File("data/pennekamp2018/K_linear-model.csv"))
df_avg = combine(
    groupby(df, [:predicted_species, :temperature, :richness, :combination]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
df_avg = innerjoin(df_avg, df_K; on = [:predicted_species, :temperature])
df_com = combine(
    groupby(df_avg, [:combination, :temperature, :richness]),
    :species_biomass => sum => :B_com,
)

inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = figure(; size = (width, 3width), fontsize = 8pt);
species_list = unique(df.predicted_species)
for (i, species) in enumerate(species_list)
    ax = Axis(fig[i, 1]; xlabel = "Temperature", ylabel = "Species biomass")
    a = subset(df_com, :combination => ByRow(==(species)))
    scatter!(a.temperature, a.community_biomass)
    i == length(species_list) || hidexdecorations!(ax)
    Label(fig[i, 2], species; rotation = 3pi / 2, tellheight = false)
end
fig

T_list = df.temperature |> unique |> sort
T_ref = first(T_list)
df_mono = subset(df_com, :richness => ByRow(==(1)))
df_poly = subset(df_com, :richness => ByRow(>(0)))
# Add reference community_biomass at temperature 15 for each combination and richness.
df_poly_ref = subset(df_poly, :temperature => ByRow(==(T_ref)))
rename!(df_poly_ref, :B_com => :B_com_ref)
select!(df_poly_ref, Not(:temperature)) # Remove temperature column for join.
df_poly = leftjoin(df_poly, df_poly_ref; on = [:combination, :richness])
# Compute difference to reference (at temperature 15)
df_poly.diff_B_com = (df_poly.B_com .- df_poly.B_com_ref) ./ df_poly.B_com_ref
df_poly = subset(df_poly, :temperature => ByRow(>(T_ref)))
df_poly.diff_B_abs = (df_poly.B_com .- df_poly.B_com_ref) ./ (df_poly.temperature .- T_ref)
df_poly.diff_B_rel = (df_poly.B_com .- df_poly.B_com_ref) ./ df_poly.B_com_ref

for a in groupby(df_poly, [:combination, :temperature, :richness])
    comb = a.combination |> first
    sp_comb = split(comb, ", ")
    sdf = subset(df_mono, :combination => ByRow(in(sp_comb)))
    B_ref = subset(sdf, :temperature => ByRow(==(T_ref))).B_com |> sum
    T_press = a.temperature |> first
    B_press = subset(sdf, :temperature => ByRow(==(T_press))).B_com |> sum
    diff_exp = (B_press - B_ref) / B_ref
    a.diff_exp = [diff_exp]
    a.s = a.diff_B_com / diff_exp
end
df_poly.diff_exp = Float64.(df_poly.diff_exp)
df_poly.normalization = df_poly.diff_B_rel ./ df_poly.s

df_poly.richness_log = log.(df_poly.richness)
df_poly.temperature_centered = df_poly.temperature .- mean(df_poly.temperature)
m = fit(
    MixedModel,
    @formula(normalization ~ temperature_centered * richness_log + (1 | combination)),
    df_poly,
)

width = 17.8cm
fig = Figure(; size = (width, 0.4width), fontsize = 10pt);
ax = Axis(fig[1, 1]; xlabel = "Richness", ylabel = "Community sensitivity")
n = nrow(df_poly)
scatter!(df_poly.richness + 0.3rand(n), df_poly.s; alpha = 0.3, color = df_poly.temperature)
hlines!([1]; color = :grey, linestyle = :dash)
ax = Axis(fig[1, 2]; xlabel = "Richness", ylabel = "Community resistance [Biomass/°C]")
scatter!(
    df_poly.richness + 0.3rand(n),
    df_poly.diff_B_abs;
    alpha = 0.5,
    color = df_poly.temperature,
)
ax = Axis(
    fig[1, 3];
    xlabel = "Pess effect on monoculture",
    ylabel = "Press effect on community",
)
temps = sort(unique(df_poly.temperature))
for t in temps
    df_t = subset(df_poly, :temperature => ByRow(==(t)))
    n = nrow(df_t)
    # Scatter points for this temperature
    scatter!(
        ax,
        df_t.diff_exp,
        df_t.diff_B_rel;
        alpha = 0.5,
        color = df_t.temperature, # vector of colors
        colorrange = (17, 25),
        label = string(t) * "°C",
    )
end
ablines!(0, 1; color = :black)
Colorbar(
    fig[1, 4];
    limits = (17, 25),
    colormap = cgrad(:viridis, 5; categorical = true),
    ticks = temps,
)
fig
