data/white2020_processed.csv: data/white2020.csv src/white2020/process.jl
	julia --project=. src/white2020/process.jl
