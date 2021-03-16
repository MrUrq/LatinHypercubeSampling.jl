module LatinHypercubeSampling

export  randomLHC,
        scaleLHC,
        AudzeEglaisObjective,
        AudzeEglaisObjective!,
        LHCoptim,
        LHCoptim!,
        subLHCoptim,
        subLHCindex,
        Categorical,
        Continuous,
        LHCDimension

using StatsBase
using Random
import Random.randperm

abstract type LHCDimension end

struct Categorical <: LHCDimension
    levels::Int64
    weight::Float64
end

struct Continuous <: LHCDimension
end

Categorical(x) = Categorical(x,0.0)

include("GA.jl")
include("AudzeEglaisObjective.jl")

#make an @threads-equivalent that is a no-op if threading is not requested
macro maybe_threaded(flag, ex)
    esc(quote
        if !$flag
            $ex
        else
            Threads.@threads $ex
        end
    end)
end

function randperm(rng,dim::Continuous,n)
    randperm(rng,n)
end

function randperm(rng,dim::Categorical,n)
    lvls = dim.levels
    out = Array{Int64}(undef,n)

    for i = 1:lvls
        s = round(Int,(i-1)*n/lvls) + 1
        e = round(Int,i*n/lvls)
        out[s:e] .= i
    end
    shuffle!(rng,out)    
end

"""
    function scaleLHC(LHC,scale_range)

Scale a Latin Hypercube to the sizes specified in scale_range. Where the range for
each dimension is specified in an vector of tuples.

# Examples
```
julia> plan = randomLHC(5,2)
5×2 Array{Int64,2}:
 5  3
 4  4
 2  2
 1  5
 3  1
julia> scaleLHC(plan,[(-1,1),(10,100)])
5×2 Array{Float64,2}:
  1.0   55.0
  0.5   77.5
 -0.5   32.5
 -1.0  100.0
  0.0   10.0
```
"""
function scaleLHC(LHC::Array{T,2},scale_range::Array{Tuple{U,V},1}) where V where U where T
    
    scaledLHC = Array{Float64,2}(undef,size(LHC))
    
    for i in 1:size(LHC,2)
        LHCcol = LHC[:,i]
        old_min, old_max = extrema(LHCcol)
        new_min, new_max = scale_range[i]

        @. scaledLHC[:,i] = (((LHCcol - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min
    end
    return scaledLHC
end

"""
    function randomLHC(n::Int, d::Int)
Generate a random Latin Hypercube with `d` dimensions and `n` sample points.
"""
function randomLHC(n::Int,d::Int)
    randomLHC(Random.GLOBAL_RNG,n,d)
end

function randomLHC(rng,n::Int,d::Int)

    dims = [Continuous() for i in 1:d]
        
    return randomLHC(rng,n,dims)

end

function randomLHC(n::Int,dims::T) where T<:AbstractArray{S} where S<:LHCDimension
    randomLHC(Random.GLOBAL_RNG,n,dims)
end

function randomLHC(rng,n::Int,dims::T) where T<:AbstractArray{S} where S<:LHCDimension

    LHC = Array{Int}(undef,n,length(dims))
    for (i, dim) = enumerate(dims)
        LHC[:,i] = randperm(rng,dim,n)
    end

    return LHC

end

"""
    function LHCoptim(n::Int,d::Int,gens;   rng::U=Random.GLOBAL_RNG,
                                            popsize::Int=100,
                                            ntour::Int=2,
                                            ptour=0.8,
                                            dims::Array{T,1}=[Continuous() for i in 1:d],
                                            interSampleWeight::Float64=1.0,
                                            periodic_ae::Bool=false,
                                            ae_power::Union{Int,Float64}=2) where T <: LHCDimension where U <: AbstractRNG
Produce an optimized Latin Hyper Cube with `d` dimensions and `n` sample points.
Optimization is run for `gens` generations. Returns a tuple of the sample plan and 
the optimization fitness history.
"""
function LHCoptim(n::Int,d::Int,gens;   rng::U=Random.GLOBAL_RNG,
                                        popsize::Int=100,
                                        ntour::Int=2,
                                        ptour=0.8,
                                        dims::Array{T,1}=[Continuous() for i in 1:d],
                                        interSampleWeight::Float64=1.0,
                                        periodic_ae::Bool=false,
                                        ae_power::Union{Int,Float64}=2,
                                        threading=false) where T <: LHCDimension where U <: AbstractRNG

    #populate first individual
    X = randomLHC(rng,n,d)

    LHCoptim!(X,gens,rng=rng,popsize=popsize,ntour=ntour,ptour=ptour,
                dims=dims,interSampleWeight=interSampleWeight,periodic_ae=periodic_ae,ae_power=ae_power,threading=threading)

end

"""
    function LHCoptim!(X::Array{Int,2},gens;    rng::U=Random.GLOBAL_RNG,
                                                popsize::Int=100,
                                                ntour::Int=2,
                                                ptour::Float64=0.8,
                                                dims::Array{T,1}=[Continuous() for i in 1:size(X,2)],
                                                interSampleWeight::Float64=1.0,
                                                periodic_ae::Bool=false,
                                                ae_power::Union{Int,Float64}=2) where T <: LHCDimension where U <: AbstractRNG
Same as LHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8) but using an
existing plan. Useful for continued optimization. Returns a tuple of the sample plan and 
the optimization fitness history.
"""
function LHCoptim!(X::Array{Int,2},gens;    rng::U=Random.GLOBAL_RNG,
                                            popsize::Int=100,
                                            ntour::Int=2,
                                            ptour::Float64=0.8,
                                            dims::Array{T,1}=[Continuous() for i in 1:size(X,2)],
                                            interSampleWeight::Float64=1.0,
                                            periodic_ae::Bool=false,
                                            ae_power::Union{Int,Float64}=2,
                                            threading=false) where T <: LHCDimension where U <: AbstractRNG

    #preallocate memory
    n, d = size(X)                             #Num points, num dimensions
    mut_inds = Array{Int}(undef,2)             #Storage of indices to swap in mutation
    tour_inds = Array{Int}(undef,ntour)        #Storage of indices for tournament selection
        

    #allocate first population
    pop = [randomLHC(rng,n,dims) for i = 1:popsize]
    pop[1] = X

    nextpop = deepcopy(pop)
    fitness = Vector{Float64}(undef, popsize)
    fitnessInds = Vector{Int64}(undef, popsize)
    bestfits = Vector{Float64}(undef, gens+1)


    #crossover for even population size
    popEven::Int64 = iseven(popsize)*-1


    #dynamic mutation rate table
    muts = Array{Int}(undef, gens)
    for l = 1:gens
        muts[l] = round(Int,-n/(gens*0.75)*l+n)
        if muts[l] < 1
            muts[l] = 1
        end
    end


    #evaluate first populations fitness
    @maybe_threaded threading for i = 1:popsize
        fitness[i] = AudzeEglaisObjective(pop[i];
                                            interSampleWeight=interSampleWeight,
                                            dims=dims,
                                            periodic_ae=periodic_ae,
                                            ae_power=ae_power)
    end


    #save the best individual and its fitness
    bestfit, bestind = findmax(fitness)
    nextpop[1] = copy(pop[bestind])
    bestfits[1] = bestfit


    #ensure fixed crossover is only applied to the continuous dimensions
    continuousDims = findall(dims.==Ref(Continuous()))

    #allocate thread-local storage for crossover
    numthreads = Threads.nthreads()
    offsprone = [similar(pop[1][:,1]) for i in 1:numthreads]
    offsprtwo = [similar(pop[1][:,1]) for i in 1:numthreads]
    randlock = ReentrantLock()

    #iterate for gens generations
    for k = 1:gens

        #tournament selection
        sortperm!(fitnessInds,fitness)
        for i = 2:popsize            
            winnerInd = tournament!(rng,fitnessInds,ntour,tour_inds,ptour)
            nextpop[i] = copy(pop[winnerInd])
        end
        
        #create children from crossover
        @maybe_threaded threading for i = 2:2:popsize+popEven
            tid = Threads.threadid()
            for j in continuousDims
                lock(randlock)
                    r = rand(rng)
                unlock(randlock)
                if r < 1.0/length(continuousDims)
                    parone = nextpop[i]
                    partwo = nextpop[i+1]
                    
                    fixedcross!(rng,offsprone[tid], offsprtwo[tid], view(parone,:,j), view(partwo,:,j),randlock)
                    nextpop[i][:,j], nextpop[i+1][:,j] = offsprone[tid], offsprtwo[tid]
                end
            end
        end

        #perform inversion and mutation
        for i = 2:popsize
            for j = 1:d
                if rand(rng) < 1.0/d     #probability for inversion
                    nextindividual = nextpop[i]
                    inversion!(rng,view(nextindividual,:,j))
                end
            end

            for j=1:muts[gens]
                mutateLHC!(rng,nextpop[i],mut_inds)
            end
        end

        #evaluate fitness
        @maybe_threaded threading for i = 1:popsize
            fitness[i] = AudzeEglaisObjective(nextpop[i];
                                                interSampleWeight=interSampleWeight,
                                                dims=dims,
                                                periodic_ae=periodic_ae,
                                                ae_power=ae_power)
        end

        #set the first individual to the best and save the fitness
        bestfit, bestind = findmax(fitness)
        nextpop[1] = nextpop[bestind]
        bestfits[k+1] = bestfit
        pop = deepcopy(nextpop)
    end

    return nextpop[1], bestfits

end

"""
    function subLHCoptim(X,n::Int,gens; rng::U=Random.GLOBAL_RNG,
                                        popsize::Int=100,
                                        ntour::Int=2,
                                        ptour::Float64=0.8,
                                        periodic_ae::Bool=false,
                                        ae_power::Union{Int,Float64}=2) where U <: AbstractRNG
Produce an optimized Latin Hyper Cube with `n` sample points from a subset of points in
`X`. Optimization is run for `gens` generations. Returns a tuple of the sample plan and
the optimization fitness history.
"""
function subLHCoptim(X,n::Int,gens; rng::U=Random.GLOBAL_RNG,
                                    popsize::Int=100,
                                    ntour::Int=2,
                                    ptour::Float64=0.8,
                                    periodic_ae::Bool=false,
                                    ae_power::Union{Int,Float64}=2) where U <: AbstractRNG

    #preallocate memory
    nLarge, d = size(X)
    pop = [Array{Int}(undef,n,d) for i = 1:popsize+1]
    nextpop = similar(pop)
    fitness = Vector{Float64}(undef, popsize+1)
    fitnessInds = Vector{Int64}(undef, popsize+1)
    bestfits = Array{Float64}(undef, gens)

    tour_inds = Array{Int}(undef,ntour)     #Storage of indices for tournament selection
    tour_fitinds = Array{Int}(undef,ntour)  #Storage of fitness for tournament selection

    #dynamic mutation rate table
    muts = Array{Int}(undef, gens)
    for l = 1:gens
        muts[l] = round(Int,-n/(gens*0.75)*l+n)
        if muts[l] < 1
            muts[l] = 1
        end
    end


    #populate first population with random subLHC's and evaluate fitness
    for i = 1:popsize+1
        subInds = sample(rng, 1:nLarge, n, replace = false)
        pop[i] = X[subInds,:]
        fitness[i] = AudzeEglaisObjective(pop[i];
                                            periodic_ae=periodic_ae,
                                            ae_power=ae_power,
                                            periodic_n=nLarge)
    end


    #save the best individual and its fitness
    bestfit, bestind = findmax(fitness)
    nextpop[1] = copy(pop[bestind])
    bestfits[1] = bestfit


    #iterate for gens generations
    for k = 1:gens
        
        #tournament selection
        sortperm!(fitnessInds,fitness)
        for i = 2:popsize+1
            winnerInd = tournament!(rng,fitnessInds,ntour,tour_inds,ptour)
            nextpop[i] = copy(pop[winnerInd])
        end
        

        #perform mutation
        for i = 2:popsize+1
            for j=1:muts[gens]
                subPoint = sample(rng,1:n) #Randomly chosen point subLHC
                largePoint = sample(rng,1:nLarge) #Randomly chosen point largeLHC

                while nextpop[i][subPoint,:] == X[largePoint,:]
                    largePoint = sample(rng,1:nLarge) #Choose new random point in largeLHC
                end
                nextpop[i][subPoint,:] = X[largePoint,:]
            end
        end

        #evaluate fitness
        for i = 1:popsize+1
            fitness[i] = AudzeEglaisObjective(nextpop[i];
                                                periodic_ae=periodic_ae,
                                                ae_power=ae_power,
                                                periodic_n=nLarge)
        end

        #set the first individual to the best and save the fitness
        bestfit, bestind = findmax(fitness)
        nextpop[1] = nextpop[bestind]
        bestfits[k] = bestfit
        pop = deepcopy(nextpop)        
    end

    return pop[1], bestfits

end

"""
    function subLHCindex(X,Xsub)
Index in the large LHC to get the subLHC.
"""
function subLHCindex(X,Xsub)

    nsub = size(Xsub,1)
    subInds = Array{Int}(undef, nsub)

    for i = 1:nsub
        A = all(Xsub[i,:]' .== X, dims=2)
        subInds[i] = (LinearIndices(A))[findall(A)][1]
    end
    return subInds
end

end # module
