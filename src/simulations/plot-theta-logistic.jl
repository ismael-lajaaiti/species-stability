using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
using LinearAlgebra
set_theme!(theme_minimal())
Random.seed!(123)

function θlogistic(B, θ, r, K)
    abs(r * (1 - (B / K)^θ))
end

r = 1
K = 2
θ = 1
B_val = LinRange(0, K, 1_000)
θ_val = LinRange(-0.1, 1, 5)

# Create the community.
S = 30
mu, sigma = -0.15, 0.1
K_std = 0.3
xi_0_avg = 0.1
D = LogNormal(log(xi_0_avg), 0.65)
df = DataFrame(; rs=Float64[], ry=Float64[], ry_hat=Float64[], r_cara=Float64[], r_short=Float64[], s=Float64[], θ=Float64[])
coms = []
n_rep = 500
for θ in θ_val
    c = rand(
        Community,
        S;
        A_ij=Normal(mu / S, sigma / sqrt(S)),
        K_i=Normal(1, K_std),
        θ=fill(θ, S),
        u=fill(1, S),
        interaction=:core,
    )
    ry = relative_yield(c)
    rs = relative_selfregulation(c)
    # Compute characteristic recovery.
    t_obs = 1
    Beq = abundance(c)
    r_cara = zeros(S)
    for i in 1:S
        x0 = zeros(S)
        x0[i] = -0.01 * Beq[i]
        sol = solve(c, Beq .+ x0, (0, t_obs))
        x_end = sol[end][i] - Beq[i]
        r_cara[i] = -1 / t_obs * log(abs.(x_end / x0[i]))
    end
    θ = c.θ
    K_hat = -c.A * Beq
    ry_hat = Beq ./ K_hat
    # Compute short-term recovery rate.
    r_short_matrix = zeros(S, n_rep)
    t_obs = 1
    for k in 1:n_rep
        x0 = -Beq .* rand(D, S)
        B0 = Beq + x0
        sol = solve(c, B0, (0, t_obs))
        x_end = sol[end] .- Beq
        r_short_matrix[:, k] = -1 / t_obs * log.(abs.(x_end ./ x0))
    end
    r_short = vec(mean(r_short_matrix, dims=2))
    # Compute sensitivity to press.
    s_matrix = zeros(S, n_rep)
    c_copy = deepcopy(c)
    for k in 1:n_rep
        du_u = -0.1 .* rand(D, S)
        c_copy.K = c.K .* (1 .+ du_u)
        Bpress = abundance(c_copy)
        s_matrix[:, k] = (Bpress .- Beq) ./ Beq ./ du_u
    end
    s = vec(mean(s_matrix, dims=2))
    append!(df, (; rs, ry, ry_hat, r_cara, r_short, s, θ))
    push!(coms, c)
end

c = coms[1]
Beq = abundance(c)
A = deepcopy(c.A)
D = Diagonal(relative_yield(c) .^ (-θ) .* Beq ./ c.K)
alpha = D * A
alpha[diagind(alpha)] .= -1
sum(inv(alpha), dims=2)

scatter((1 ./ df.ry .- a .* (1 ./ df.ry .- 1)), df.s, color=df.θ)

z0 = rand(D, 10_000)
a = mean(z0) / harmmean(z0)
scatter(df.rs .* (1 .+ a .* (1 ./ df.ry_hat .- 1)) .* df.θ, df.r_short, color=df.θ)

inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 17.8cm
fig = Figure(; size=(width, width), fontsize=10pt);
ga = fig[1, :] = GridLayout()
gb = fig[2, :] = GridLayout()
gc = fig[3, :] = GridLayout()
ax = Axis(ga[1, 2], xlabel="Biomass (B)", ylabel="Per capita growth rate (f)")
for θ in θ_val
    f_val = θlogistic.(B_val, θ, r, K)
    lines!(B_val, f_val, label="θ = $θ")
end
axislegend(; fontsize=10pt)
ax2 = Axis(gb[1, 1]; xlabel="Relative yield (RY)", ylabel="Self-regulation loss (SL)")
ry_val = LinRange(0, 1.2, 1_000)
for θ in θ_val
    dfθ = subset(df, :θ => ByRow(==(θ)), :ry => ByRow(>(0.01)))
    scatter!(dfθ.ry, dfθ.rs, label="θ=$θ")
    lines!(ry_val, ry_val .^ θ)
end
ax3 = Axis(gb[1, 2]; xlabel="Relative yield (RY)", ylabel="Characteristic recovery rate")
ax4 = Axis(gb[1, 3]; xlabel="Self-regulation loss (SL)", ylabel="Characteristic recovery rate")
ax5 = Axis(gc[1, 1]; xlabel="Relative yield (RY)", ylabel="Pseudo relative yield (RY)")
ax6 = Axis(gc[1, 2]; xlabel="Short-term  recovery rate", ylabel="Prediction based on RY")
ax7 = Axis(gc[1, 3]; xlabel="Sensitivity to press", ylabel="Prediction based on RY")
for θ in θ_val
    @info θ
    dfθ = subset(df, :θ => ByRow(==(θ)), :ry => ByRow(>(0.01)))
    scatter!(ax3, dfθ.ry, dfθ.r_cara, label="θ=$θ")
    scatter!(ax4, dfθ.rs, dfθ.r_cara, label="θ=$θ")
    scatter!(ax5, dfθ.ry, dfθ.ry_hat)
    r_short_pred = dfθ.rs .* (1 .+ a .* (1 ./ dfθ.ry .- 1))
    scatter!(ax6, dfθ.r_short, r_short_pred)
    s_pred = 1 ./ dfθ.ry .- a .* (1 ./ dfθ.ry .- 1)
    scatter!(ax7, dfθ.s, s_pred)
end
ablines!(ax3, 0, 1, color=:black)
ablines!(ax4, 0, 1, color=:black)
ablines!(ax5, 0, 1, color=:black)
ablines!(ax6, 0, 1, color=:black)
ablines!(ax7, 0, 1, color=:black)
# hideydecorations!(ax4)
fig

save("figures/simulations/thetalogistic.svg", fig)
