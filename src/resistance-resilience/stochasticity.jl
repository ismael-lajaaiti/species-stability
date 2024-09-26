import DifferentialEquations: DiscreteCallback, CallbackSet
using GLV
using DataFrames
using Distributions
using Random
using CairoMakie
include("makie-theme.jl")
include("scripts/resistance-resilience/functions.jl")
set_theme!(theme_minimal())
Random.seed!(1234)

S = 30
mu, sigma = -1, 0.3
K_std = 0.3
c = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
)
N_eq = abundance(c)

noise_intensity = 0.06
noise!(du, u, p, t) =
    for i in eachindex(du)
        du[i] = noise_intensity * u[i]
    end

# Press parameters.
K_no_press = deepcopy(c.K)
K_press = deepcopy(K_no_press)
beta = 0.1 # Press intensity.
i = 1 # Index of the species undergoing the press.
K_press[i] = (1 - beta) * K_no_press[i]
press_on!(integrator) = integrator.p.K = K_press
press_off!(integrator) = integrator.p.K = K_no_press
t_start, t_stop = 100, 200
press_start(u, t, integrator) = t in [t_start]
press_stop(u, t, integrator) = t in [t_stop]
cb_start = DiscreteCallback(press_start, press_on!)
cb_stop = DiscreteCallback(press_stop, press_off!)
callback = CallbackSet(cb_start, cb_stop)
tstops = [t_start, t_stop]
sol = solve(c, abundance(c), (1, 300), noise!; callback, tstops)

fig = Figure();
ax = Axis(fig[1, 1]; xlabel = "Time", ylabel = "Abundance")
lines!(sol.t, sol[1, :]; color = :black, linewidth = 0.5)
fig

df1, _ = get_resistance_resilience(
    c;
    noise_intensity,
    stochastic = true,
    beta = 0.1,
    n_experiments = 1,
)
df10, _ = get_resistance_resilience(
    c;
    noise_intensity,
    stochastic = true,
    beta = 0.1,
    n_experiments = 10,
)

size = (700, 300)
fig = Figure(; size);
ax1 = Axis(
    fig[1, 1];
    xlabel = "Resistance",
    ylabel = "Resilience",
    xscale = log10,
    yscale = log10,
    title = "1 experiment",
)
scatter!(df1.resistance, df1.resilience; color = df1.relative_yield)
ax2 = Axis(
    fig[1, 2];
    xlabel = "Resistance",
    xscale = log10,
    yscale = log10,
    title = "10 experiments",
)
scatter!(df10.resistance, df10.resilience; color = df10.relative_yield)
limits = extrema(df1.relative_yield)
linkaxes!(ax1, ax2)
fig[1, 3] = Colorbar(fig[1, 2]; label = "Relative yield", limits)
fig
save_figure("figures/resist-resil-experiments", fig, size)
