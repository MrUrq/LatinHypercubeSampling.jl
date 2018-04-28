# LatinHypercubeSampling.jl

| **Documentation** | **Build & Testing Status** |
|:-----------------:|:--------------------------:|
[![][docs-stable-img]][docs-stable-url] | [![Build Status](https://travis-ci.org/MrUrq/LatinHypercubeSampling.jl.svg?branch=master)](https://travis-ci.org/MrUrq/LatinHypercubeSampling.jl) [![Coverage Status](https://coveralls.io/repos/github/MrUrq/LatinHypercubeSampling.jl/badge.svg?branch=master)](https://coveralls.io/github/MrUrq/LatinHypercubeSampling.jl?branch=master) [![codecov.io](http://codecov.io/github/MrUrq/LatinHypercubeSampling.jl/coverage.svg?branch=master)](http://codecov.io/github/MrUrq/LatinHypercubeSampling.jl?branch=master) | 

*LatinHypercubeSampling* is a Julia package for the creation of optimised Latin Hypercube Sampling Plans. The genetic optimisation algorithim is largely based on the work by Bates et al. [1]. The package includes additional functionality for the creation of an optimised subset of an existing plan.

Features:

* Creation of an optimised Latin Hypercube Sampling plan.
* Generate an optimised subset of an existing plan.
* Refine existing plan through mutation only.




[1]: Stuart Bates, Johann Sienz, and Vassili Toropov. "Formulation of the Optimal Latin Hypercube Design of Experiments Using a Permutation Genetic Algorithm", 45th AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics & Materials Conference, Structures, Structural Dynamics, and Materials and Co-located Conferences, () https://doi.org/10.2514/6.2004-2011

Picture of LHC

## Examples

A quite extensive set of examples can be found at the [PGFPlotsXExamples repo](https://github.com/KristofferC/PGFPlotsXExamples).

## Installation

The package is registered in `METADATA.jl` and can be installed with `Pkg.add`.

```julia
julia> Pkg.add("LatinHypercubeSampling")
```

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **tagged version of the documentation.**


## Author

- Magnus Urquhart - [@MrUrq](https://github.com/MrUrq/)

[docs-stable-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://MrUrq.github.io/LatinHypercubeSampling.jl/latest
