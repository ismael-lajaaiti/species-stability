# Stability of species within ecological communities

## Abstract

> [!CAUTION]
> Add abstract.

## Reproduce figures of the article

To reproduce figures of the article, first you have to clone this repository.

Secondly, you can move to the `scripts/` directory and install the required package automatically with

```bash
cd scripts # From the root of the project.
julia --project=. -e 'using Pkg; Pkg.instantiate()' # Install Julia packages.
```

> [!TIP]
> It is advised to use Julia 1.10, otherwise you may need to resolve conflicts between package versions with
> ```bash
> julia --project=. -e 'using Pkg; Pkg.update()'
> ```

Finally, you can execute the scripts with

```bash
julia --project=. <name_of_the_script>.jl
```

where `<name_of_the_script>` is the name of the script you want to run.
Figures are saved in the `scripts/figures/` directory.


## Structure of the code

This project is structured as a Julia package named `LotkaVolterra`.
The package contains general functions used by the different scripts producing the figures of the article.
The code of the package is stored in the `LotkaVolterra/` directory.
The scripts producing the figures of the article are stored in the `scripts/` directory.
The figures produced by the scripts are stored in `scripts/figures/` directory (which is created after the first script is run).
