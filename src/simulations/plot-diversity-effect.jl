using GLV
using CairoMakie
using LinearAlgebra
using DataFrames
using Distributions

# Create the community.
S = 30
sigma = 0.1
mu = 0.1
K_std = 0.3
c = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
)
B = abundance(c)
K = c.K
rho_val = LinRange(-0.3, 1, 10)
kappa_avg = 0.1
n_rep = 100

df = DataFrame(; rho = Float64[], s = Float64[])
for rho in rho_val
    for k in 1:n_rep
        c_copy = deepcopy(c)
        kappa = 0.1 ./ c.K #generate_kappa_with_cov(c.K, rho) .* kappa_avg
        c_copy.K = c.K .* (1 .+ kappa)
        deltaK = c_copy.K - c.K
        Bpress = abundance(c_copy)
        deltaB = Bpress .- B
        s = sum(deltaB) / sum(B) * sum(K) / sum(deltaK)
        push!(df, (rho, s))
    end
end
combine(groupby(df, :rho), :s => mean)

V = -inv(c.A)
Vii = V[diagind(V)] |> mean
Vij = sum(V - Diagonal(V)) / (S * (S - 1))
rho = 1.2
s = 1 / (1 + rho) + rho / (1 + rho) * sum(K) / sum(B) * (Vii + (S - 1) * Vij)

sensitivity(deltaB, deltaK, B, K) = (deltaB / B) * (K / deltaK)
df = select(df, :, [:deltaB, :deltaK, :B, :K] => ByRow(sensitivity) => :s)
df_com =
    combine(groupby(df, [:rho, :p]), [:B, :deltaB, :K, :deltaK] .=> sum; renamecols = false)
df_com =
    select(df_com, :, [:deltaB, :B] => ByRow(/) => :dB, [:deltaK, :K] => ByRow(/) => :dK)
df_com = select(df_com, :rho, [:dB, :dK] => ByRow(/) => :s_com)
df_com = combine(groupby(df_com, :rho), :s_com => mean; renamecols = false)


# Create the community.
S = 30
sigma = 0.1
mu_val = LinRange(-0.15, 0.1, 20)
K_std = 0.3
xi_0_avg = 0.1
rho = -2.0
kappa_avg = 0.01
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
coms = []
t_end = 1 # Observation time to compute recovery rates.
n_rep = 500
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
        kappa = generate_kappa_with_cov(c.K, rho) .* kappa_avg
        c_copy.K = c.K .* (1 .+ kappa)
        deltaK = c_copy.K - c.K
        Bpress = abundance(c_copy)
        deltaB = Bpress .- B
        p = fill(k, S)
        append!(df, (; p, mu, sp, deltaB, deltaK, B, K))
    end
    push!(coms, c)
end

function generate_kappa_with_cov(x::Vector{Float64}, target_cov::Float64; noise_std = 0.1)
    n = length(x)
    a = target_cov / var(x)
    z = randn(n) .* noise_std
    a .* x .+ z
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

function expected_dB(kappa_mean, K_mean, cov_kappa_K, mu, S)

end
