using Statistics
using CSV
using DataFrames
using GLM
using CairoMakie
using StatsBase
using Distributions
set_theme!(theme_minimal())

# Process data.
day_start = 7 # Let some time for community to stabilise.
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

# Plot.
fig = Figure(; size = (600, 550));
ax1 = Axis(fig[1, 1]; xlabel = "Carrying capacity (μg/mL)", ylabel = "Biomass (μg/mL)")
ax2 = Axis(fig[1, 2]; xlabel = "Scaled carrying capacity", ylabel = "Scaled biomass")
ax3 = Axis(
    fig[2, 1:2];
    xlabel = "Species relative yield",
    ylabel = "Species sensitivity",
    aspect = AxisAspect(1.5),
)
level = 0.9 # Confidence interval for linear model.
df_s =
    DataFrame(; ry = Float64[], s_mean = Float64[], e_low = Float64[], e_high = Float64[])
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
    scatter!(
        ax1,
        gdf.K_estimated,
        gdf.species_biomass;
        label = "$sp",
        markersize = gdf.temperature .- 8,
        color = B_ref / K_ref,
        colormap = :algae,
        colorrange = (0, 1),
        alpha = 0.5,
    )
    lines!(
        ax1,
        gdf.K_estimated,
        gdf.predicted_biomass;
        color = B_ref / K_ref,
        colormap = :algae,
        colorrange = (0, 1),
    )
    scatter!(
        ax2,
        gdf.K_estimated ./ K_ref,
        gdf.species_biomass ./ B_ref;
        label = "$sp",
        alpha = 0.5,
        color = B_ref / K_ref,
        colormap = :algae,
        colorrange = (0, 1),
        markersize = gdf.temperature .- 8,
    )
    lines!(
        ax2,
        gdf.K_estimated ./ K_ref,
        gdf.predicted_biomass ./ B_ref;
        color = B_ref / K_ref,
        colormap = :algae,
        colorrange = (0, 1),
    )
end
for row in eachrow(df_s)
    scatter!(ax3, row.ry, row.s_mean)
    errorbars!(ax3, [row.ry], [row.s_mean], [row.e_low], [row.e_high]; whiskerwidth = 10)
end
ry_min, ry_max = extrema(df_s.ry)
ry = LinRange(ry_min, ry_max, 100)
s_ii = 1 ./ ry
pred_sensitivity = (s_ii .+ (1 .- s_ii) * am_on_hm_avg)
lines!(ax3, ry, pred_sensitivity; color = :black, label = "analytical prediction")
axislegend()
markersizes = unique(df_avg.temperature) .- 8
group_size = [
    MarkerElement(;
        marker = :circle,
        color = :black,
        strokecolor = :transparent,
        markersize = ms,
    ) for ms in markersizes
]
fig[1, 3] = Legend(fig, group_size, string.(markersizes .+ 8), "Temperature")
fig[2, 3] = Legend(fig, ax1, "Species")
l1 = fig[1, 1] = GridLayout()
l2 = fig[1, 2] = GridLayout()
l3 = fig[2, 1:2] = GridLayout()
for (label, layout) in zip(["A", "B", "C"], [l1, l2, l3])
    Label(
        layout[1, 1, TopLeft()],
        label;
        font = :bold,
        padding = label == "C" ? (0, -20, 5, 0) : (0, 5, 5, 0),
        halign = :right,
    )
end
fig

save("figures/pennekamp2018/resistance-vs-ry.png", fig)
save("figures/pennekamp2018/resistance-vs-ry.svg", fig)
