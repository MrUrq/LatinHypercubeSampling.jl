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
Where `gens` is the number of generations to run the optimisation for. The number of
samples selected for tournament in the genetic optimisation algorithm as well as the
probability for tournament selection can be accessed with the optional arguments `ntour=2`
and `ptour=0.8`. These are manually tuned defaults and may not be optimal depending on the
number of dimensions and the size of the plan.

The optimisation of a sampling plan is started from a random plan which is also
an exported function.
```@docs
randomLHC(n::Int, d::Int)
```

## Example
The `LHCoptim` function run for many generations to create an optimised 120 point
plan in 2 dimensions. Both the plan and the fitness of the sampling plan are exported.

```julia-repl
julia> plan, fitness = LHCoptim(120,2,gens)
```

If only the plan is required it can be omitted using 
```julia-repl
julia> plan, _ = LHCoptim(120,2,gens)
```

![](https://raw.githubusercontent.com/MrUrq/LatinHypercubeSampling.jl/master/docs/src/assets/120p2d.png)
