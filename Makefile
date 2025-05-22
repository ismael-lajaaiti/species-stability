.PHONY: setup all main supporting # High-level actions.

all: setup main supporting # Build everything.

setup:
	mkdir -p figures
	julia --project=. -e 'using Pkg; Pkg.instantiate()'

main: data/pennekamp2018/processed-data.csv \
      data/pennekamp2018/K_linear-model.csv \
	  data/pennekamp2018/A_normalized.csv \
      figures/data.png \
	  figures/pulse.png \
	  figures/press.png

supporting: figures/si-theta-logistic.png \
			figures/si-competition-gradient.png \
			figures/si-feedback.png \
			figures/si-dependent-species.png \
			figures/si-sensitivity-mu.png

data/pennekamp2018/processed-data.csv: data/pennekamp2018/raw-data.csv src/pennekamp2018/process.jl
	julia --project=. src/pennekamp2018/process.jl

data/pennekamp2018/K_linear-model.csv: data/pennekamp2018/processed-data.csv src/pennekamp2018/compute-K.jl
	julia --project=. src/pennekamp2018/plot-biomass-vs-richness.jl

data/pennekamp2018/A_normalized.csv: data/pennekamp2018/processed-data.csv src/pennekamp2018/infer-interactions.jl
	julia --project=. src/pennekamp2018/infer-interactions.jl

figures/data.png: data/pennekamp2018/processed-data.csv data/pennekamp2018/K_linear-model.csv src/pennekamp2018/plot-data.jl
	julia --project=. src/pennekamp2018/plot-data.jl

figures/pulse.png: src/simulations/plot-pulse.jl
	julia --project=. src/simulations/plot-pulse.jl

figures/press.png: src/simulations/plot-press.jl
	julia --project=. src/simulations/plot-press.jl

figures/si-theta-logistic.png: src/simulations/plot-theta-logistic.jl
	julia --project=. src/simulations/plot-theta-logistic.jl

figures/si-competition-gradient.png: src/simulations/plot-competition-gradient.jl
	julia --project=. src/simulations/plot-competition-gradient.jl

figures/si-feedback.png: src/simulations/plot-feedback.jl
	julia --project=. src/simulations/plot-feedback.jl

figures/si-dependent-species.png: src/simulations/plot-dependent-species.jl
	julia --project=. src/simulations/plot-dependent-species.jl

figures/si-sensitivity-mu.png: src/simulations/plot-sensitivity-mu.jl
	julia --project=. src/simulations/plot-sensitivity-mu.jl
