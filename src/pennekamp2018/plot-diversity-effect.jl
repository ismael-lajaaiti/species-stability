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
day_start = 1 # Let some time for community to stabilise.
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
df_A = DataFrame(CSV.File("data/pennekamp2018/A_normalized.csv"))

T_list = df.temperature |> unique |> sort
T_ref = first(T_list)
df_mono = subset(df_com, :richness => ByRow(==(1)))
df_poly = subset(df_com, :richness => ByRow(>=(1)))
# Add reference community_biomass at temperature 15 for each combination and richness.
df_poly_ref = subset(df_poly, :temperature => ByRow(==(T_ref)))
rename!(df_poly_ref, :B_com => :B_com_ref)
select!(df_poly_ref, Not(:temperature)) # Remove temperature column for join.
df_poly = leftjoin(df_poly, df_poly_ref; on = [:combination, :richness])
# Compute difference to reference (at temperature 15)
df_poly.diff_B_com = (df_poly.B_com .- df_poly.B_com_ref) ./ df_poly.B_com_ref
df_poly = subset(df_poly, :temperature => ByRow(>(T_ref)))
df_poly.diff_B_abs = (df_poly.B_com .- df_poly.B_com_ref) ./ (df_poly.temperature .- T_ref)
df_poly.diff_B_rel =
    (df_poly.B_com .- df_poly.B_com_ref) ./ df_poly.B_com_ref ./
    (df_poly.temperature .- T_ref)

for a in groupby(df_poly, [:combination, :temperature, :richness])
    comb = a.combination |> first
    sp_comb = split(comb, ", ")
    sdf = subset(df_mono, :combination => ByRow(in(sp_comb)))
    B_ref = subset(sdf, :temperature => ByRow(==(T_ref))).B_com |> sum
    B_ref_sp = subset(sdf, :temperature => ByRow(==(T_ref))).B_com
    T_press = a.temperature |> first
    B_press = subset(sdf, :temperature => ByRow(==(T_press))).B_com |> sum
    diff_exp = (B_press - B_ref) / B_ref #/ (T_press - T_ref)
    a.diff_exp = [diff_exp]
    a.s = a.diff_B_com ./ diff_exp
end
df_poly.diff_exp = Float64.(df_poly.diff_exp)
# df_poly.normalization = df_poly.diff_B_rel ./ df_poly.s

for a in groupby(df_poly, [:combination, :temperature, :richness])
    S = a.richness |> first
    if S != 1
        t = a.temperature |> first
        df_A_t = subset(df_A, :temperature => ByRow(==(t)))
        comb = a.combination |> first
        sp_comb = split(comb, ", ")
        subset!(df_A_t, :row => ByRow(in(sp_comb)))
        select!(df_A_t, sp_comb)
        a.a = [sum(Array(df_A_t)) / (S * (S - 1))]
    else
        a.a = [0]
    end
end
df_poly.a = Float64.(df_poly.a)

df_sp = combine(
    groupby(df_avg, [:combination, :temperature, :richness, :predicted_species]),
    :species_biomass => :B_sp,
)
df_sp_ref = subset(df_sp, :temperature => ByRow(==(T_ref)))
rename!(df_sp_ref, :B_sp => :B_ref)
select!(df_sp_ref, Not(:temperature))
df_sp = leftjoin(df_sp, df_sp_ref; on = [:combination, :predicted_species, :richness])
df_sp = subset(df_sp, :temperature => ByRow(>(T_ref)))
df_sp.diff_B_rel = (df_sp.B_sp .- df_sp.B_ref) ./ df_sp.B_ref
df_mono_diff = subset(df_mono, :temperature => ByRow(>(T_ref)))
df_mono_ref = subset(df_mono, :temperature => ByRow(==(T_ref)))
select!(df_mono_ref, Not(:temperature))
rename!(df_mono_ref, :B_com => :B_ref)
df_mono_diff = leftjoin(df_mono_diff, df_mono_ref; on = [:combination, :richness])
df_mono_diff.exp = (df_mono_diff.B_com .- df_mono_diff.B_ref) ./ df_mono_diff.B_ref
rename!(df_mono_diff, :combination => :predicted_species)
select!(df_mono_diff, Not(:richness, :B_com, :B_ref))
df_sp = innerjoin(df_sp, df_mono_diff; on = [:predicted_species, :temperature])
df_sp.s = df_sp.diff_B_rel ./ df_sp.exp

for a in groupby(df_sp, [:combination, :temperature, :richness, :predicted_species])
    S = a.richness |> first
    if S != 1
        t = a.temperature |> first
        sp = a.predicted_species |> first
        df_A_t = subset(df_A, :temperature => ByRow(==(t)))
        comb = a.combination |> first
        sp_comb = split(comb, ", ")
        subset!(df_A_t, :row => ByRow(==(sp)))
        select!(df_A_t, sp_comb)
        @info df_A_t
        a.a = [sum(Array(df_A_t)) / (S - 1)]
    else
        a.a = [0]
    end
end
df_sp.a = Float64.(df_sp.a)

# Compute species coefficient of variation (CV) for each combination and temperature
df_cv = combine(
    groupby(df, [:combination, :temperature, :predicted_species]),
    :species_biomass => (x -> std(skipmissing(x)) / mean(skipmissing(x))) => :cv,
)
subset!(df_cv, :temperature => ByRow(>(15)))
df_sp = innerjoin(df_sp, df_cv; on = [:temperature, :combination, :predicted_species])

scatter(1 ./ df_sp.cv, df_sp.diff_B_rel; color = df_sp.richness)

# Fit a mixed effects model to test the relationship between cv and diff_B_rel
m_cv = fit(MixedModel, @formula(diff_B_rel ~ cv + (1 | combination)), df_sp)
m_cv

# Community scale
# Compute community biomass per replicate, then community-level CV for each combination and temperature
df_com_biomass = combine(
    groupby(df, [:combination, :temperature, :day, :replicate]),
    :species_biomass => sum => :community_biomass,
)
df_cv_com = combine(
    groupby(df_com_biomass, [:combination, :temperature, :replicate]),
    :community_biomass =>
        (x -> mean(skipmissing(x)) / std(skipmissing(x))) => :invariability,
)
df_cv_com = combine(
    groupby(df_cv_com, [:combination, :temperature]),
    :invariability => mean => :invariability,
)
subset!(df_cv_com, :temperature => ByRow(>(15)))
df_poly = innerjoin(df_poly, df_cv_com; on = [:temperature, :combination])

scatter(df_poly.invariability, df_poly.diff_B_abs; color = df_poly.richness)
m_cv = fit(MixedModel, @formula(diff_B_abs ~ invariability + (1 | combination)), df_poly)
m_cv

# df_poly.diff_B_abs = abs.(df_poly.diff_B_abs)
# df_poly.diff_B_rel = abs.(df_poly.diff_B_rel)
# df_poly.s = abs.(df_poly.s)
df_poly.richness_log = log.(df_poly.richness)
df_poly.temperature_centered = df_poly.temperature .- mean(df_poly.temperature)
m1 = fit(
    MixedModel,
    @formula(diff_B_abs ~ temperature_centered * richness + (1 | combination)),
    df_poly,
)
@info m1.pvalues
m2 = fit(
    MixedModel,
    @formula(diff_B_rel ~ temperature_centered * richness + (1 | combination)),
    df_poly,
)
m3 = fit(MixedModel, @formula(diff_B_rel ~ a * temperature + (1 | combination)), df_poly)

df_sp.s_abs = abs.(df_sp.s)
df_sp.richness_log = log.(df_sp.richness)
df_sp.temperature_centered = df_sp.temperature .- mean(df_sp.temperature)
m = fit(
    MixedModel,
    @formula(s_abs ~ temperature_centered * richness_log + (1 | combination)),
    df_sp,
)


inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
width = 17.8cm
alpha = 0.8
fig = Figure(; size = (width, 0.4width), fontsize = 10pt);
ax = Axis(
    fig[1, 1];
    xlabel = "Richness",
    ylabel = "Absolute community\nbiomass change[μg/mL/°C]",
)
n = nrow(df_poly)
m = nrow(df_sp)
scatter!(
    df_poly.richness + 0.3rand(n),
    df_poly.diff_B_abs;
    alpha,
    color = df_poly.temperature,
)
if m1.pvalues[3] < 0.05
    a = m1.pvalues[1] < 0.05 ? coef(m1)[1] : 0
    b = coef(m1)[3]
    ablines!(a, b; color = :black, linewidth = 2)
end
ax = Axis(
    fig[1, 2];
    xlabel = "Richness",
    ylabel = "Relative community\nbiomass change [/°C]",
    # yscale = log10,
)
scatter!(df_sp.a, abs.(df_sp.s); alpha, color = df_sp.temperature)
ax = Axis(
    fig[1, 3];
    xlabel = "Richness",
    ylabel = "1 + biotic contribution to\nbiomass change",
    yscale = log10,
)
scatter!(ax, df_poly.a, abs.(df_poly.s); alpha, color = df_poly.richness)
# ylims!(ax, [-1.5, 10])
# temps = sort(unique(df_poly.temperature))
# for t in temps
#     df_t = subset(df_poly, :temperature => ByRow(==(t)))
#     n = nrow(df_t)
#     # Scatter points for this temperature
#     scatter!(
#         ax,
#         df_t.richness .+ 0.3rand(n),
#         df_t.s;
#         alpha,
#         color = df_t.temperature, # vector of colors
#         colorrange = (17, 25),
#         label = string(t) * "°C",
#     )
# end
# hlines!([0]; color = :grey, linestyle = :dash, linewidth = 2)
Colorbar(
    fig[1, 4];
    limits = extrema(df_poly.a),
    # colormap = cgrad(:viridis, 5; categorical = true),
    # ticks = temps,
)
fig
