# LatinHypercubeSampling

## Introduction
*LatinHypercubeSampling* is a Julia package for the creation of optimised Latin Hypercube (LHC) Sampling Plans. The genetic optimisation algorithm is largely based on the work by Bates et al. [1]. The package includes additional functionality for the creation of an optimised subset of an existing plan.

## Installation
The package is registered and can be installed with `Pkg.add`.

```julia
julia> Pkg.add("LatinHypercubeSampling")
```

### Reference
[1]: Stuart Bates, Johann Sienz, and Vassili Toropov. "Formulation of the Optimal Latin Hypercube Design of Experiments Using a Permutation Genetic Algorithm", 45th AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics & Materials Conference, Structures, Structural Dynamics, and Materials and Co-located Conferences, () https://doi.org/10.2514/6.2004-2011
