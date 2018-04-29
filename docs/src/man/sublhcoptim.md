# Optimised subset of LHC Sampling Plan

Generate an optimised subset of an existing plan. The optimisation of the subset is
based on a genetic algorithm.

## Functions
```@docs
subLHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)
```

```@docs
subLHCindex(X,Xsub)
```

## Example
Create an optimised subset LHC plan from an existing plan with 120 points in 2 dimensions.

![](https://raw.githubusercontent.com/MrUrq/LatinHypercubeSampling.jl/master/docs/src/assets/120p2d.png)

```julia-repl
julia> subLHCoptim(X,Xsub)
```

![](https://raw.githubusercontent.com/MrUrq/LatinHypercubeSampling.jl/master/docs/src/assets/sub60p2d.png)

The indices of the subset in the larger plan can be extracted with
```julia-repl
julia> subLHCindex(X,Xsub)
```
