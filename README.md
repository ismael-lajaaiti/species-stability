# Revealing the organization of species stability in ecological communities

Code to reproduce the figures of the manuscript *Revealing the organization of species stability in ecological communities*.

## Generate the figures of the article

To reproduce figures of the article, you first have to clone this repository.
Then move to the cloned repository and run in a terminal
The Julia environment can be setup with

```sh
make setup
```

This should create the `figures/` folder, and install all required packages.

> [!TIP]
> We advise to use Julia 1.11.

Next, the figures of the main text can be generate with

```sh
make main
```

Finally, to generate the figures of the supporting information, run

```sh
make supporting
```
