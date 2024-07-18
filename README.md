## paleoAM

`paleoAM` is an R package for fitting models of abundance to individual species from fossil assemblages, and simulating those assemblages under various ecological and geological biases. It includes functions for modeling abundance distributions over environmental gradients, using kernal density estimation, as well as tools for simulating communities that might be sampled at different positions along a gradient, and how our recover of the community would vary under different sedimentary histories and sampling approaches.

Users of `paleoAM` functions can fit models of abundance to their data, and then simulate how time-averaging, bioturbation and other mixing processes might change or obscure a given paleoenvironmental signal that they infer from their ecological data.
	
You can install this latest development version using the R function `install_github` in the package `devtools`:

```
devtools::install_github("dwbapst/paleoAM")
```
	
Once installed, you can check the version number of your `paleotree` install using the R function `packageVersion`:

```
packageVersion("paleoAM")
```

The following badge indicates the package Status of R CHECK via GitHub actions:

[![R Check](https://github.com/dwbapst/paleoAM/actions/workflows/r.yml/badge.svg)](https://github.com/dwbapst/paleoAM/actions)
  <!-- badges: end -->

