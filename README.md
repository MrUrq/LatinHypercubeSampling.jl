<img src="docs/src/assets/logo.png" width="180">

# LatinHypercubeSampling.jl

| **Documentation** | **Build & Testing Status** |
|:-----------------:|:--------------------------:|
[![][docs-stable-img]][docs-stable-url] | [![build](https://github.com/MrUrq/LatinHypercubeSampling.jl/workflows/CI/badge.svg)](https://github.com/MrUrq/LatinHypercubeSampling.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/MrUrq/LatinHypercubeSampling.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MrUrq/LatinHypercubeSampling.jl) | 

*LatinHypercubeSampling* is a Julia package for the creation of optimised Latin Hypercube Sampling Plans. The genetic optimisation algorithm is largely based on the work by Bates et al. [1]. The package includes additional functionality for the creation of an optimised subset of an existing plan. For more details, see our [paper](https://doi.org/10.1016/j.asoc.2019.106050).

Features:

* Creation of an optimised Latin Hypercube Sampling plan.
* Generate an optimised subset of an existing plan.
* Refine existing plan.
* Ability to include discrete parameters in the design.

It also has the option to optimize the sampling plans using the periodic Audze–Eglājs criteria [2].

## Installation

The package is registered and can be installed with `Pkg.add`.

```julia
julia> Pkg.add("LatinHypercubeSampling")
```

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **tagged version of the documentation.**


## Author

- Magnus Urquhart - [@MrUrq](https://github.com/MrUrq/)

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://MrUrq.github.io/LatinHypercubeSampling.jl/stable

## Example 
Sampling the Rosenbrock function with an optimized Latin Hypercube sampling plan.
```julia-repl
julia> plan, _ = LHCoptim(100,2,1000)
julia> scaled_plan = scaleLHC(plan,[(-5.0,5.0),(-5.0,5.0)])
julia> rosenbrock_2D(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
julia> mapslices(rosenbrock_2D,scaled_plan; dims=2)
```

## Example LHC
Example of optimised LHC plan for 120 points in 2 dimensions.
<img src="docs/src/assets/120p2d.png">

### References
[1]: Stuart Bates, Johann Sienz, and Vassili Toropov. "Formulation of the Optimal Latin Hypercube Design of Experiments Using a Permutation Genetic Algorithm", 45th AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics & Materials Conference, Structures, Structural Dynamics, and Materials and Co-located Conferences, () https://doi.org/10.2514/6.2004-2011

[2]: Jan Eliáš, Miroslav Vořechovský, Modification of the Audze–Eglājs criterion to achieve a uniform distribution of sampling points, Advances in Engineering Software, Volume 100, 2016, Pages 82-96, ISSN 0965-9978, () https://doi.org/10.1016/j.advengsoft.2016.07.004.

### Citation
```
@article{urquhart_surrogate-based_2020,
	title = {Surrogate-based optimisation using adaptively scaled radial basis functions},
	volume = {88},
	issn = {1568-4946},
	doi = {10.1016/j.asoc.2019.106050},
	journal = {Applied Soft Computing},
	author = {Urquhart, Magnus and Ljungskog, Emil and Sebben, Simone},
	year = {2020},
}
```
