using GLV
using StatsBase
using DataFrames
using Distributions
using Random
using CairoMakie
using LinearAlgebra
set_theme!(theme_minimal())
Random.seed!(123)


# Create the community.
S = 30
sigma = 0.1
mu_val = LinRange(-0.15, 0, 20)
K_std = 0.3
xi_0_avg = 0.1
D = LogNormal(log(xi_0_avg), 0.65)
df = DataFrame(;
    p = Float64[],
    mu = Float64[],
    sp = Int64[],
    deltaB = Float64[],
    deltaK = Float64[],
    B = Float64[],
    K = Float64[],
)
df2 =
    DataFrame(; p = Float64[], mu = Float64[], sp = Int64[], xt = Float64[], x0 = Float64[])
coms = []
t_end = 1 # Observation time to compute recovery rates.
n_rep = 100
for mu in mu_val
    c = rand(
        Community,
        S;
        A_ij = Normal(mu / S, sigma / sqrt(S)),
        K_i = Normal(1, K_std),
        interaction = :core,
    )
    B = abundance(c)
    K = c.K
    sp = collect(1:S)
    # Compute sensitivity to press.
    c_copy = deepcopy(c)
    mu = fill(mu, S)
    for k in 1:n_rep
        du_u = -0.1 .* rand(D, S)
        c_copy.K = c.K .* (1 .+ du_u)
        deltaK = c_copy.K - c.K
        Bpress = abundance(c_copy)
        deltaB = Bpress .- B
        p = fill(k, S)
        append!(df, (; p, mu, sp, deltaB, deltaK, B, K))
    end
    # Short-term recovery to pulse.
    for k in 1:n_rep
        x0 = -0.1 * rand(D, S) .* B
        sol = solve(c, B .+ x0, (0, t_end))
        xt = sol[end] .- B
        p = fill(k, S)
        append!(df2, (; p, mu, sp, xt, x0))
    end
    push!(coms, c)
end

sensitivity(deltaB, deltaK, B, K) = (deltaB / B) * (K / deltaK)
df = select(df, :, [:deltaB, :deltaK, :B, :K] => ByRow(sensitivity) => :s)
df_sp = combine(groupby(df, :mu), :s => mean => :s_sp)
df_com =
    combine(groupby(df, [:mu, :p]), [:B, :deltaB, :K, :deltaK] .=> sum; renamecols = false)
df_com =
    select(df_com, :, [:deltaB, :B] => ByRow(/) => :dB, [:deltaK, :K] => ByRow(/) => :dK)
df_com = select(df_com, :mu, [:dB, :dK] => ByRow(/) => :s_com)
df_com = combine(groupby(df_com, :mu), :s_com => mean; renamecols = false)
df_join = innerjoin(df_sp, df_com; on = :mu)

recovery_rate(xt, x0) = -log(abs(xt / x0)) / t_end
df2 = select(df2, :, [:xt, :x0] => ByRow(recovery_rate) => :r)
df_sp2 = combine(groupby(df2, :mu), :r => mean => :r_sp)
df_com2 = combine(groupby(df2, [:mu, :p]), :xt => sum, :x0 => sum)
df_com2 = select(df_com2, :, [:xt_sum, :x0_sum] => ByRow(recovery_rate) => :r_com)
df_com2 = combine(groupby(df_com2, :mu), :r_com => mean => :r_com)
df_join2 = innerjoin(df_sp2, df_com2; on = :mu)

inch = 96
pt = 4 / 3
cm = inch / 2.54
width = 11cm
fig = Figure(; size = (width, 0.6 * width), fontsize = 10pt);
ax = Axis(
    fig[1, 1];
    xlabel = "Mean interaction strength",
    ylabel = "Sensitivity to press \n (reversed)",
)
ax.yreversed = true
scatter!(df_join.mu, df_join.s_com; label = "community", color = :grey)
scatter!(df_join.mu, df_join.s_sp; label = "species", color = :black)
ax = Axis(
    fig[1, 2];
    xlabel = "Mean interaction strength",
    ylabel = "Short-term recover rate",
)
scatter!(df_join2.mu, df_join2.r_com; label = "community", color = :grey)
scatter!(df_join2.mu, df_join2.r_sp; label = "species", color = :black)
fig[0, :] = axislegend(; orientation = :horizontal, halign = :center)
fig

save("figures/si-competition-gradient.png", fig)
# save("figures/simulations/competition-gradient.svg", fig)
