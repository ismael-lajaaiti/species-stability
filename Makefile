.PHONY: setup all clean main supporting # High-level actions.

all: setup main supporting # Build everything.

clean:
	rm -rf figures
	rm -f data/pennekamp2018/processed-data.csv
	rm -f data/pennekamp2018/A_normalized.csv
	rm -f data/pennekamp2018/df_mono.csv

setup:
	mkdir -p figures
	julia --project=. -e 'using Pkg; Pkg.instantiate()'

main: data/pennekamp2018/processed-data.csv \
	  data/pennekamp2018/A_normalized.csv \
	  data/pennekamp2018/df_mono.csv \
      figures/data.png \
	  figures/pulse.png \
	  figures/press.png

supporting: figures/si-theta-logistic.png \
			figures/si-competition-gradient.png \
			figures/si-feedback.png \
			figures/si-dependent-species.png \
			figures/si-sensitivity-mu.png \
			figures/si-strong-interactions.png \
			figures/si-press-others.png \
			figures/si-eigvec-alignment.png \
			figures/si-pulse-growth-rate.png

data/pennekamp2018/processed-data.csv: data/pennekamp2018/raw-data.csv src/pennekamp2018/process.jl
	julia --project=. src/pennekamp2018/process.jl

data/pennekamp2018/A_normalized.csv data/pennekamp2018/df_mono.csv: data/pennekamp2018/processed-data.csv src/pennekamp2018/infer-A-and-K.jl
	julia --project=. src/pennekamp2018/infer-A-and-K.jl

figures/data.png: data/pennekamp2018/processed-data.csv data/pennekamp2018/df_mono.csv data/pennekamp2018/A_normalized.csv src/pennekamp2018/plot-data.jl
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

figures/si-strong-interactions.png: src/simulations/plot-strong-interactions.jl
	julia --project=. src/simulations/plot-strong-interactions.jl

figures/si-press-others.png: src/simulations/plot-press-others.jl
	julia --project=. src/simulations/plot-press-others.jl

figures/si-eigvec-alignment.png: src/simulations/plot-eigvec-alignment.jl
	julia --project=. src/simulations/plot-eigvec-alignment.jl

figures/si-pulse-growth-rate.png: src/simulations/plot-growth-rate.jl
	julia --project=. src/simulations/plot-growth-rate.jl