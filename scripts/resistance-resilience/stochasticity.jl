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

noise_intensity = 0.05
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
