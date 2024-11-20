using Statistics
using CSV
using DataFrames
using GLM
using CairoMakie
using StatsBase
set_theme!(theme_minimal())

# Analysis parameters.
day_start = 7 # Let some time for community to stabilise.
S_min = 3 # Minimal commuunity richness, for analysis.

# Process data.
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
df = subset(df, :day => ByRow(>=(day_start)))
data = combine(
    groupby(df, [:predicted_species, :richness, :temperature]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
df_K = DataFrame(CSV.File("data/pennekamp2018/K_linear-model.csv"))
df_all = subset(df, :richness => ByRow(>=(S_min)))
df_all = combine(
    groupby(df_all, [:predicted_species, :temperature]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
df_all = innerjoin(df_all, df_K; on = [:predicted_species, :temperature])

df_kappa = DataFrame(; predicted_species = String[], kappa = Float64[])
for gdf in groupby(df_all, :predicted_species)
    sp = first(unique(gdf.predicted_species))
    Kmin, Kmax = extrema(gdf.K_estimated)
    kappa = (Kmax - Kmin) / mean(gdf.K_estimated)
    push!(df_kappa, (sp, kappa))
end

kappa_min, kappa_max = extrema(df_kappa.kappa)
# mu = mean(rand(Uniform(kappa_min, kappa_max), 1_000_000))
mu_reciprocal = mean(1 ./ rand(Uniform(kappa_min, kappa_max), 1_000_000))
mu = mean(df_kappa.kappa)
# mu_reciprocal = mean(1 ./ df_kappa.kappa)

ry_min, ry_max = extrema(df_stab.ry)
ry = LinRange(ry_min, ry_max, 100)
s_ii = 1 ./ ry
pred_stab = 1 ./ (s_ii .+ (1 .- s_ii) * mu * mu_reciprocal)

# Plot.
fig = Figure(; size = (500, 500));
ax1 = Axis(fig[1, 1]; xlabel = "K", ylabel = "B")
ax2 = Axis(fig[1, 2]; xlabel = "K / <K>", ylabel = "B / <B>")
ax3 = Axis(fig[2, 1:2]; xlabel = "η", ylabel = "Species stability \n to press")
df_stab = DataFrame(; ry = Float64[], stab = Float64[])
for gdf in groupby(df_all, :predicted_species)
    sp = gdf.predicted_species |> first
    B_ref = mean(gdf.species_biomass)
    K_ref = mean(gdf.K_estimated)
    model = lm(@formula(species_biomass ~ K_estimated), gdf)
    s = coef(model)[2]
    stab = abs((B_ref / K_ref) * (1 / s))
    push!(df_stab, (B_ref / K_ref, stab))
    gdf.predicted_biomass = predict(model, gdf)
    gdf = dropmissing(gdf, :predicted_biomass)
    scatter!(ax1, gdf.K_estimated, gdf.species_biomass; label = "$sp", alpha = 0.8)
    lines!(ax1, gdf.K_estimated, gdf.predicted_biomass)
    scatter!(
        ax2,
        gdf.K_estimated ./ K_ref,
        gdf.species_biomass ./ B_ref;
        label = "$sp",
        alpha = 0.8,
    )
    lines!(ax2, gdf.K_estimated ./ K_ref, gdf.predicted_biomass ./ B_ref)
end
scatter!(ax3, df_stab.ry, df_stab.stab; color = :black)
lines!(ax3, ry, pred_stab)
fig[0, 1:2] = Legend(fig, ax1, "Species"; orientation = :horizontal)
fig

save("figures/pennekamp2018/resistance-vs-ry.png", fig)
