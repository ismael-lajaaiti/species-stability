.PHONY: all main supporting # High-level actions.

all: main supporting # Build everything.

main: data/pennekamp2018/processed-data.csv \
      data/pennekamp2018/K_linear-model.csv \
      figures/data.png \
	  figures/pulse.png \
	  figures/press.png

data/pennekamp2018/processed-data.csv: data/pennekamp2018/raw-data.csv src/pennekamp2018/process.jl
	julia --project=. src/pennekamp2018/process.jl

data/pennekamp2018/K_linear-model.csv: data/pennekamp2018/processed-data.csv src/pennekamp2018/plot-biomass-vs-richness.jl
	julia --project=. src/pennekamp2018/plot-biomass-vs-richness.jl

figures/data.png: data/pennekamp2018/processed-data.csv data/pennekamp2018/K_linear-model.csv src/pennekamp2018/plot-resistance-vs-ry.jl
	julia --project=. src/pennekamp2018/plot-resistance-vs-ry.jl

figures/pulse.png: src/simulations/plot-pulse.jl
	julia --project=. src/simulations/plot-pulse.jl

figures/press.png: src/simulations/plot-press.jl
	julia --project=. src/simulations/plot-press.jl
