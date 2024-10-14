# df_K = DataFrame(CSV.File("data/pennekamp2018/carrying-capacity.csv"))
df = subset(df, :richness => ByRow(==(6)))
df = subset(df, :microcosmID => ByRow(==(697)))
T = df.temperature |> unique |> first
df_K = subset(df_K, :temperature => ByRow(==(T)))

data = DataFrame(; species = String[], variance = Float64[], B = Float64[])
for g in groupby(df, :predicted_species)
    species = g.predicted_species |> first
    v = var(g.species_biomass)
    b = mean(g.species_biomass)
    push!(data, (species, v, b))
end
data = innerjoin(data, df_K; on = :species)
data.eta = data.B ./ data.K

fig = Figure();
ax = Axis(fig[1, 1]; xlabel = "Biomass", ylabel = "CV", yscale = log10, xscale = log10)
scatter!(data.B, sqrt.(data.variance))
fig
