using Statistics
using CSV
using DataFrames
using GLM
using CairoMakie
using StatsBase
using Distributions
set_theme!(theme_minimal())

# Process data.
day_start = 10 # Let some time for community to stabilise.
S_min = 3 # Minimal community richness, for analysis.
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
df = subset(df, :day => ByRow(>=(day_start)))
df = subset(df, :richness => ByRow(>=(S_min)))
df_K = DataFrame(CSV.File("data/pennekamp2018/K_linear-model.csv"))
df_avg = combine(
    groupby(df, [:predicted_species, :temperature]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
df_avg = innerjoin(df_avg, df_K; on = [:predicted_species, :temperature])

df_avg2 = combine(
    groupby(df, [:predicted_species, :temperature, :combination]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
df_avg2 = innerjoin(df_avg2, df_K; on = [:predicted_species, :temperature])

# Compute arithmetic and harmonic of press perturbation intensities.
df_kappa = DataFrame(; predicted_species = String[], kappa = Float64[])
for gdf in groupby(df_avg, :predicted_species)
    sp = first(unique(gdf.predicted_species))
    Kmin, Kmax = extrema(gdf.K_estimated)
    kappa = (Kmax - Kmin) / mean(gdf.K_estimated)
    push!(df_kappa, (sp, kappa))
end
kappa_list = df_kappa.kappa
am_on_hm = zeros(6)
for i in eachindex(kappa_list)
    am_on_hm[i] = harmmean(vcat(kappa_list[begin:i-1], kappa_list[i+1:end])) / kappa_list[i]
end
am_on_hm_avg = mean(am_on_hm)

level = 0.8 # Confidence interval for linear model.
df_s2 = DataFrame(; ry = Float64[], s_mean = Float64[], e = [], species = [])
for gdf in groupby(df_avg2, [:predicted_species, :combination])
    B_ref = mean(gdf.species_biomass)
    K_ref = mean(gdf.K_estimated)
    ry_ref = B_ref / K_ref
    ry = mean(gdf.species_biomass ./ gdf.K_estimated)
    model = lm(@formula(species_biomass ~ K_estimated), gdf)
    s_mean = coef(model)[2] / ry_ref
    S = unique(gdf.combination) |> first
    s_low = coeftable(model; level).cols[5][2] / ry_ref
    s_high = coeftable(model; level).cols[6][2] / ry_ref
    e = s_mean - s_low
    species = gdf.predicted_species |> unique |> first
    push!(df_s2, (ry, s_mean, e, species))
end
df_s2

# Plot.
inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, width), fontsize = 8pt);
l1 = fig[1, 1:3] = GridLayout()
l2 = fig[2, 1:3] = GridLayout()
ax1 = Axis(l2[1, 1]; xlabel = "Carrying capacity (μg/mL)", ylabel = "Biomass (μg/mL)")
ax2 = Axis(l2[1, 2]; xlabel = "Carrying capacity (μg/mL)")
hideydecorations!(ax2)
ax3 = Axis(
    l1[1, 1];
    xlabel = "SL",
    ylabel = "Sensitivity to press\n(reversed)",
    # aspect=AxisAspect(1.5),
)
df_s =
    DataFrame(; ry = Float64[], s_mean = Float64[], e_low = Float64[], e_high = Float64[])
colorrange = extrema(df.temperature)
colormap = :lipari
for gdf in groupby(df_avg, :predicted_species)
    sp = gdf.predicted_species |> first
    B_ref = mean(gdf.species_biomass)
    K_ref = mean(gdf.K_estimated)
    ry_ref = B_ref / K_ref
    ry = mean(gdf.species_biomass ./ gdf.K_estimated)
    model = lm(@formula(species_biomass ~ K_estimated), gdf)
    s_mean = coef(model)[2] / ry_ref
    s_low = coeftable(model; level).cols[5][2] / ry_ref
    s_high = coeftable(model; level).cols[6][2] / ry_ref
    e_low = s_mean - s_low
    e_high = s_high - s_mean
    push!(df_s, (ry, s_mean, e_low, e_high))
    gdf.predicted_biomass = predict(model, gdf)
    gdf = dropmissing(gdf, :predicted_biomass)
    scatter!(ax1, gdf.K_estimated, gdf.species_biomass; label = "$sp", alpha = 0.5)
    lines!(ax1, gdf.K_estimated, gdf.predicted_biomass;)
    lines!(ax2, gdf.K_estimated, gdf.predicted_biomass; color = :grey)
    scatter!(
        ax2,
        gdf.K_estimated,
        gdf.species_biomass;
        label = "$sp",
        alpha = 0.7,
        color = gdf.temperature,
        colorrange,
        colormap,
    )
end

colors = Makie.wong_colors()
species = df.predicted_species |> unique
color_dict = Dict(sp => color for (color, sp) in zip(colors, species))
for sp in species
    color = color_dict[sp]
    df_sp = subset(df_s2, :species => ByRow(==(sp)))
    scatter!(
        ax3,
        df_sp.ry,
        df_sp.s_mean;
        markersize = 12 .- 8 .* df_sp.e,
        alpha = 0.7,
        color = color_dict[sp],
    )
end
ry_min, ry_max = extrema(df_s2.ry)
ry = LinRange(ry_min + 0.05, ry_max, 100)
s_ii = 1 ./ ry
pred_sensitivity = (s_ii .+ (1 .- s_ii) * am_on_hm_avg)
lines!(ax3, ry, pred_sensitivity; color = :black, label = "analytical prediction")
elems = [LineElement(), MarkerElement(; marker = :circle)]
axislegend(ax3, elems, ["analytical\nprediction", "data"]; position = :rt)
ax3.yreversed = true
l1[1, 2] = Legend(fig, ax1, "Species"; rowgap = -4, tellwidth = true)
cb = Colorbar(l2[1, 3]; limits = colorrange, colormap, label = "Temperature (°C)")
for (label, layout) in zip(["A", "B", "C"], [l1, l2[1, 1], l2[1, 2]])
    Label(
        layout[1, 1, TopLeft()],
        label;
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right,
    )
end
fig

save("figures/data.png", fig)
# save("figures/simulations/data.svg", fig)
# save("figures/pennekamp2018/resistance-vs-ry.png", fig)
# save("figures/pennekamp2018/resistance-vs-ry.svg", fig)
