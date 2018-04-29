# Optimised Latin Hypercube Sampling Plan

Create an optimised Latin Hypercube Sampling Plan using a genetic based optimisation
algorithm. The objective function is the inverse of the Audze-Eglais function defined as

```math
\text{max } U = \text{max} \sum_{p=1}^P \sum_{q=p+1}^P L^2_{pq}
```
where ``L^2_{pq}`` is the square of the Euclidean norm.
!!! note

    This package maximises the inverse of the Audze-Eglais objective function.

## Function
```@docs
LHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)
```

## Example
The `LHCoptim` function run for many generations to create an optimised 120 point
plan in 2 dimensions.

```julia-repl
julia> LHCoptim(120,2,gens)
```

![](https://raw.githubusercontent.com/MrUrq/LatinHypercubeSampling.jl/master/docs/src/assets/120p2d.png)
