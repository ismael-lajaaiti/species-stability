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
day_end = 250
S_min = 3 # Minimal community richness, for analysis.
df = DataFrame(CSV.File("data/pennekamp2018/processed-data.csv"))
df = subset(df, :day => ByRow(>=(day_start)))
df = subset(df, :day => ByRow(<=(day_end)))
df = subset(df, :richness => ByRow(>=(S_min)))
df_avg = combine(
    groupby(df, [:predicted_species, :temperature]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
df_A = DataFrame(CSV.File("data/pennekamp2018/A_normalized.csv"))
df_mono = DataFrame(CSV.File("data/pennekamp2018/df_mono.csv"))
df_avg = innerjoin(df_avg, df_mono; on = [:predicted_species, :temperature])
rename!(df_avg, :B_mono => :K)
df_avg2 = combine(
    groupby(df, [:predicted_species, :temperature, :combination]),
    :species_biomass => mean ∘ skipmissing => :species_biomass,
)
df_avg2 = innerjoin(df_avg2, df_mono; on = [:predicted_species, :temperature])

df_SL = deepcopy(df_avg2)
subset!(df)
rename!(df_SL, :species_biomass => :B, :B_mono => :K)
df_SL.SL = df_SL.B ./ df_SL.K
gdf = groupby(df_SL, [:temperature, :combination])
for a in gdf
    sp_list = a.predicted_species
    SL_exp_list = []
    Ai_list = []
    T = a.temperature |> first
    df_A_T = subset(df_A, :temperature => ByRow(==(T)))
    for row in eachrow(a)
        focal_sp = row.predicted_species
        A_incoming = subset(df_A_T, :row => ByRow(==(focal_sp)))
        Ai = 0
        SL_exp = 1
        mean_SL = 0
        for interacting_sp in sp_list
            if interacting_sp != focal_sp
                SL_j = a[a.predicted_species.==interacting_sp, :SL] |> first
                mean_SL += SL_j / (length(sp_list) - 1)
                mean_SL += SL_exp += A_incoming[1, interacting_sp] * SL_j
                Ai += A_incoming[1, interacting_sp]
            end
        end
        push!(SL_exp_list, SL_exp)
        push!(Ai_list, Ai)
    end
    a.SL_exp = SL_exp_list
    a.Ai = Ai_list
end

species_list = unique(df.predicted_species)
df_mu = DataFrame(; species = String[], temperature = Float64[], mu = Float64[])
for gdf in groupby(df_A, :temperature)
    T = gdf.temperature |> first
    gdf_sp = select(gdf, species_list)
    mu0 = [sum(col) / 5 for col in eachcol(gdf_sp)]
    for (i, species) in enumerate(names(gdf_sp))
        push!(df_mu, (species, T, mu0[i]))
    end
end
df_mu = combine(groupby(df_mu, :species), :mu => mean => :mu)

# Compute arithmetic and harmonic of press perturbation intensities.
df_kappa = DataFrame(; predicted_species = String[], kappa = Float64[])
for gdf in groupby(df_avg, :predicted_species)
    sp = first(unique(gdf.predicted_species))
    Kmin, Kmax = extrema(gdf.K)
    kappa = (Kmax - Kmin) / mean(gdf.K)
    push!(df_kappa, (sp, kappa))
end
kappa_list = df_kappa.kappa
am_on_hm = zeros(6)
for i in eachindex(kappa_list)
    am_on_hm[i] = harmmean(vcat(kappa_list[begin:i-1], kappa_list[i+1:end])) / kappa_list[i]
end
am_on_hm_avg = mean(am_on_hm)

level = 0.8 # Confidence interval for linear model.
df_s2 = DataFrame(;
    ry = Float64[],
    s_mean = Float64[],
    e = [],
    species = [],
    combination = String[],
)
for gdf in groupby(df_avg2, [:predicted_species, :combination])
    species = gdf.predicted_species |> unique |> first
    combination = gdf.combination |> unique |> first
    B_ref = mean(gdf.species_biomass)
    K_ref = mean(gdf.B_mono)
    ry_ref = B_ref / K_ref
    ry = mean(gdf.species_biomass ./ gdf.B_mono)
    mu0 = df_mu[df_mu.species.==species, :mu][1]
    Vii = (1 + mu0 * ry) / (1 + mu0)
    model = lm(@formula(species_biomass ~ B_mono), gdf)
    s_mean = coef(model)[2] / ry_ref
    S = unique(gdf.combination) |> first
    s_low = coeftable(model; level).cols[5][2] / ry_ref
    s_high = coeftable(model; level).cols[6][2] / ry_ref
    e = s_mean - s_low
    push!(df_s2, (ry, s_mean, e, species, combination))
end
df_s2

# Plot 1 - Main text.
alpha = 0.8
inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 10cm
fig = Figure(; size = (width, 1.4width), fontsize = 8pt);
l1 = fig[1, 1:3] = GridLayout()
l2 = fig[2, 1:3] = GridLayout()
ax1 = Axis(l2[1, 1]; xlabel = "Carrying capacity (μg/mL)", ylabel = "Biomass (μg/mL)")
ax2 = Axis(l2[1, 2]; xlabel = "Carrying capacity (μg/mL)")
hideydecorations!(ax2)
ax3 = Axis(l1[1, 1]; xlabel = "SL", ylabel = "Sensitivity to press\n(reversed)")
df_s =
    DataFrame(; ry = Float64[], s_mean = Float64[], e_low = Float64[], e_high = Float64[])
colorrange = extrema(df.temperature)
colormap = :lipari
for gdf in groupby(df_avg, :predicted_species)
    sp = gdf.predicted_species |> first
    B_ref = mean(gdf.species_biomass)
    K_ref = mean(gdf.K)
    ry_ref = B_ref / K_ref
    ry = mean(gdf.species_biomass ./ gdf.K)
    model = lm(@formula(species_biomass ~ K), gdf)
    s_mean = coef(model)[2] / ry_ref
    s_low = coeftable(model; level).cols[5][2] / ry_ref
    s_high = coeftable(model; level).cols[6][2] / ry_ref
    e_low = s_mean - s_low
    e_high = s_high - s_mean
    push!(df_s, (ry, s_mean, e_low, e_high))
    gdf.predicted_biomass = predict(model, gdf)
    gdf = dropmissing(gdf, :predicted_biomass)
    scatter!(ax1, gdf.K, gdf.species_biomass; label = "$sp", alpha = 0.5)
    lines!(ax1, gdf.K, gdf.predicted_biomass;)
    lines!(ax2, gdf.K, gdf.predicted_biomass; color = :grey)
    scatter!(
        ax2,
        gdf.K,
        gdf.species_biomass;
        label = "$sp",
        alpha,
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
        abs.(df_sp.s_mean);
        markersize = 12 .- 8 .* df_sp.e,
        alpha,
        color = color_dict[sp],
    )
end
ry_min, ry_max = extrema(df_s2.ry)
ry = LinRange(ry_min + 0.01, ry_max, 100)
s_ii = 1 ./ ry
pred_sensitivity = (s_ii .+ (1 .- s_ii) * am_on_hm_avg)
lines!(ax3, ry, pred_sensitivity; color = :black, label = "analytical prediction")
hlines!(1; color = :grey, linestyle = :dash, label = "baseline")
elems = [LineElement(), MarkerElement(; marker = :circle)]
axislegend(ax3, elems, ["analytical\nprediction", "data"]; position = :rt)
ax3.yreversed = true
l1[1, 2] = Legend(fig, ax1, "Species"; rowgap = -4, tellwidth = true)
cb = Colorbar(l2[1, 3]; limits = colorrange, colormap, label = "Temperature (°C)")
# Last panel - interactions.
l3 = fig[3, :] = GridLayout()
ax3 = Axis(l3[1, 1]; xlabel = "Incoming direct interactions", ylabel = "SL observed")
ax4 = Axis(l3[1, 2]; xlabel = "Incoming net interactions", ylabel = "SL observed")
for (i, sp) in enumerate(unique(df.predicted_species))
    df_sp = subset(df_SL, :predicted_species => ByRow(==(sp)))
    scatter!(ax3, df_sp.Ai, df_sp.SL; alpha)
    scatter!(ax4, 1 .+ df_sp.SL_exp, df_sp.SL; alpha)
end
hideydecorations!(ax4)
for (label, layout) in
    zip(["A", "B", "C", "D", "E"], [l1, l2[1, 1], l2[1, 2], l3[1, 1], l3[1, 2]])
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

# Prepare plot - PhD Thesis.
df_sp_a = combine(groupby(df_SL, [:predicted_species, :combination]), :SL_exp => mean)
rename!(df_sp_a, :predicted_species => :species)
df_sp = innerjoin(df_s2, df_sp_a; on = [:species, :combination])
df_com = combine(
    groupby(df, [:combination, :richness, :temperature, :day]),
    :species_biomass => sum => :B,
)
df_com = combine(groupby(df_com, [:combination, :richness, :temperature]), :B => mean => :B)
B_mono_list = []
A_list = []
for row in eachrow(df_com)
    sp_list = split(row.combination, ", ")
    T = row.temperature
    df_tmp = subset(
        df_mono,
        :temperature => ByRow(==(T)),
        :predicted_species => ByRow(in(sp_list)),
    )
    push!(B_mono_list, sum(df_tmp.B_mono))
    df_tmp = subset(df_A, :row => ByRow(in(sp_list)), :temperature => ByRow(==(T)))
    select!(df_tmp, sp_list)
    S = length(sp_list)
    A = sum(Array(df_tmp)) / S / (S - 1)
    push!(A_list, A)
end
df_com.B_mono = B_mono_list
df_com.A = A_list
df_com.B_mono = Float64.(df_com.B_mono)
for gdf in groupby(df_com, :combination)
    model = lm(@formula(B ~ B_mono), gdf)
    s = coef(model)[2] / mean(gdf.B ./ gdf.B_mono)
    gdf.s = fill(s, nrow(gdf))
end
df_com = combine(groupby(df_com, [:combination]), :s => mean, :A => mean)
df_com2 = combine(groupby(df_sp, :combination), :ry => mean)
df_com = innerjoin(df_com, df_com2; on = :combination)

alpha = 0.8
inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 13cm
fig = Figure(; size = (width, 0.45width), fontsize = 10pt);
ax = Axis(
    fig[1, 1];
    xlabel = "Mean interaction",
    ylabel = "Interaction contribution\nto resistance",
    title = "Species",
)
scatter!(df_sp.SL_exp_mean .- 1, 1 .- df_sp.s_mean; color = :black)
ax2 = Axis(fig[1, 2]; xlabel = "Mean interaction", title = "Community")
scatter!(df_com.ry_mean .- 1, 1 .- df_com.s_mean; color = :black)
fig

save("figures/fading.png", fig)

# Check relationship significance.
model = lm(@formula(s_mean ~ ry_mean), df_com)
model = lm(@formula(s_mean ~ SL_exp_mean), df_sp)
m_com = fit(MixedModel, @formula(s_mean ~ A_mean + (1 | combination)), df_com)
