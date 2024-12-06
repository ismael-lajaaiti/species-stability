.PHONY: all pennekamp2018 white2020 # High-level actions.

all: simulations pennekamp2018 white2020 # Build everything.

pennekamp2018: data/pennekamp2018/processed-data.csv \
               data/pennekamp2018/carrying-capacity.csv \
               figures/pennekamp2018/monocultures.png \
			   figures/pennekamp2018/carrying-capacity.png \
               figures/pennekamp2018/var-vs-biomass.png \
               figures/pennekamp2018/fit-logistic.png \
               figures/pennekamp2018/biomass-vs-richness.png \
               figures/pennekamp2018/resistance-vs-ry.png

simulations: figures/simulations/press-species-community.png \
	         figures/simulations/return-time.png

data/pennekamp2018/processed-data.csv: data/pennekamp2018/raw-data.csv src/pennekamp2018/process.jl
	julia --project=. src/pennekamp2018/process.jl

data/pennekamp2018/carrying-capacity.csv: data/pennekamp2018/processed-data.csv src/pennekamp2018/compute-carrying-capacity.jl
	julia --project=. src/pennekamp2018/compute-carrying-capacity.jl

data/pennekamp2018/K_linear-model.csv: data/pennekamp2018/processed-data.csv src/pennekamp2018/plot-biomass-vs-richness.jl
	julia --project=. src/pennekamp2018/plot-biomass-vs-richness.jl

figures/pennekamp2018/%.png: data/pennekamp2018/processed-data.csv src/pennekamp2018/plot-%.jl
	julia --project=. src/pennekamp2018/plot-$*.jl

figures/simulations/%.png: src/simulations/plot-%.jl
	julia --project=. src/simulations/plot-$*.jl

white2020: data/white2020_processed.csv

data/white2020_processed.csv: data/white2020.csv src/white2020/process.jl
	julia --project=. src/white2020/process.jl
