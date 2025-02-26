using GLV
using DataFrames
using Distributions
using Random
using CairoMakie
include("makie-theme.jl")
set_theme!(theme_minimal())
Random.seed!(1234)

S = 30
mu, sigma = 0.0, 0.3
K_std = 0.3
c = rand(
    Community,
    S;
    A_ij=Normal(mu / S, sigma / sqrt(S)),
    K_i=Normal(1, K_std),
    interaction=:core,
)

function get_noise!(du, u, p, t; exponent, noise_intensity)
    for i in eachindex(u)
        du[i] = noise_intensity * u[i]^exponent
    end
end

N_eq = abundance(c)
tspan = (0, 1_000)
noise_intensity = 0.05
scaling_dict = Dict("Immigration" => 0, "Demographic" => 0.5, "Environmental" => 1)
df = DataFrame(; noise_type=String[], species_id=Float64[], species_cv=Float64[])
for (noise_type, exponent) in scaling_dict
    noise!(du, u, p, t) = get_noise!(du, u, p, t; exponent, noise_intensity)
    sol = solve(c, N_eq, tspan, noise!)
    species_cv = get_species_cv(sol)
    species_id = 1:S
    df_tmp = DataFrame(; noise_type=fill(noise_type, S), species_id, species_cv)
    append!(df, df_tmp)
end

# DataFrame storing species information.
df_species = DataFrame(;
    species_id=1:S,
    abundance=abundance(c),
    relative_yield=relative_yield(c),
)
df_immigration = df[df.noise_type.=="Immigration", :]
df_immigration = innerjoin(df_immigration, df_species; on=:species_id)
df_demographic = df[df.noise_type.=="Demographic", :]
df_demographic = innerjoin(df_demographic, df_species; on=:species_id)
df_environmental = df[df.noise_type.=="Environmental", :]
df_environmental = innerjoin(df_environmental, df_species; on=:species_id)

function cv_prediction(noise_intensity, exponent, abundance, relative_yield)
    noise_intensity / sqrt(2) * abundance^(exponent - 1) / sqrt(relative_yield)
end

size = (400, 300)
alpha = 1
fig = Figure(; size);
# CV vs. relative yield.
ax =
    Axis(fig[1, 1]; xlabel="CV prediction", ylabel="CV", xscale=log10, yscale=log10)
ablines!(ax, 0, 1; color=:black)
scatter!(
    cv_prediction.(
        noise_intensity,
        0,
        df_immigration.abundance,
        df_immigration.relative_yield,
    ),
    df_immigration.species_cv;
    label="Immigration",
    alpha,
)
scatter!(
    cv_prediction.(
        noise_intensity,
        0.5,
        df_demographic.abundance,
        df_demographic.relative_yield,
    ),
    df_demographic.species_cv;
    label="Demographic",
    alpha,
)
scatter!(
    cv_prediction.(
        noise_intensity,
        1,
        df_environmental.abundance,
        df_environmental.relative_yield,
    ),
    df_environmental.species_cv;
    label="Environmental",
    alpha,
)
axislegend(ax; position=:lt)
fig
save_figure("figures/cv-prediction", fig, size)
