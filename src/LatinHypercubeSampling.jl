module LatinHypercubeSampling

export  randomLHC,
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


function randperm(dim::Continuous,n)
    randperm(n)
end

function randperm(dim::Categorical,n)
    lvls = dim.levels
    out = Array{Int64}(undef,n)

    for i = 1:lvls
        s = round(Int,(i-1)*n/lvls) + 1
        e = round(Int,i*n/lvls)
        out[s:e] .= i
    end
    shuffle!(out)    
end


"""
    function randomLHC(n::Int, d::Int)
Generate a random Latin Hypercube with `d` dimensions and `n` sample points.
"""
function randomLHC(n::Int,d::Int)

    dims = [Continuous() for i in 1:d]
        
    return randomLHC(n,dims)

end


function randomLHC(n::Int,dims::T) where T<:AbstractArray{S} where S<:LHCDimension

    LHC = Array{Int}(undef,n,length(dims))
    for (i, dim) = enumerate(dims)
        LHC[:,i] = randperm(dim,n)
    end

    return LHC

end




"""
    function LHCoptim(n::Int,d::Int,gens)
Produce an optimized Latin Hyper Cube with `d` dimensions and `n` sample points.
Optimization is run for `gens` generations.
"""
function LHCoptim(n::Int,d::Int,gens;   popsize::Int=100,
                                        ntour::Int=2,
                                        ptour=0.8,
                                        dims::Array{T,1}=[Continuous() for i in 1:d],
                                        interSampleWeight::Float64=1.0) where T <: LHCDimension

    #populate first individual
    X = randomLHC(n,d)

    LHCoptim!(X,gens,popsize=popsize,ntour=ntour,ptour=ptour,
                dims=dims,interSampleWeight=interSampleWeight)

end



"""
    function LHCoptim!(X::Array{Int,2},n::Int,d::Int,gens;ntour::Int=2,ptour=0.8)
Same as LHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8) but using an
existing population. Useful for continued optimization.
"""
function LHCoptim!(X::Array{Int,2},gens;    popsize::Int=100,
                                            ntour::Int=2,
                                            ptour::Float64=0.8,
                                            dims::Array{T,1}=[Continuous() for i in 1:size(X,2)],
                                            interSampleWeight::Float64=1.0) where T <: LHCDimension

    #preallocate memory
    n, d = size(X)                             #Num points, num dimensions
    mut_inds = Array{Int}(undef,2)             #Storage of indices to swap in mutation
    tour_inds = Array{Int}(undef,ntour)        #Storage of indices for tournament selection
        

    #allocate first population
    pop = [randomLHC(n,dims) for i = 1:popsize]
    pop[1] = X

    nextpop = deepcopy(pop)
    fitness = Vector{Float64}(undef, popsize)
    fitnessInds = Vector{Int64}(undef, popsize)
    offsprone = similar(pop[1][:,1])
    offsprtwo = similar(pop[1][:,1])
    bestfits = Vector{Float64}(undef, gens+1)
    dist = zeros(Float64,Int(n*(n-1)*0.5))


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
    for i = 1:popsize
        fitness[i] = AudzeEglaisObjective!(dist,  pop[i];
                                            interSampleWeight=interSampleWeight,
                                            dims=dims)
    end


    #save the best individual and its fitness
    bestfit, bestind = findmax(fitness)
    nextpop[1] = copy(pop[bestind])
    bestfits[1] = bestfit


    #ensure fixed crossover is only applied to the continuous dimensions
    continuousDims = findall(dims.==Ref(Continuous()))

    #iterate for gens generations
    for k = 1:gens

        #tournament selection
        sortperm!(fitnessInds,fitness)
        for i = 2:popsize            
            winnerInd = tournament!(fitnessInds,ntour,tour_inds,ptour)
            nextpop[i] = copy(pop[winnerInd])
        end
        
        #create children from crossover
        for i = 2:2:popsize+popEven
            for j in continuousDims
                if rand() < 1.0/length(continuousDims)
                    parone = nextpop[i]
                    partwo = nextpop[i+1]
                    
                    fixedcross!(offsprone, offsprtwo, view(parone,:,j), view(partwo,:,j))
                    nextpop[i][:,j], nextpop[i+1][:,j] = offsprone, offsprtwo
                end
            end
        end

        #perform inversion and mutation
        for i = 2:popsize
            for j = 1:d
                if rand() < 1.0/d     #probability for inversion
                    nextindividual = nextpop[i]
                    inversion!(view(nextindividual,:,j))
                end
            end

            for j=1:muts[gens]
                mutateLHC!(nextpop[i],mut_inds)
            end
        end

        #evaluate fitness
        for i = 1:popsize
            fitness[i] = AudzeEglaisObjective!(dist, nextpop[i];
                                                interSampleWeight=interSampleWeight,
                                                dims=dims)
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
    function subLHCoptim(X,n::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)
Produce an optimized Latin Hyper Cube with `n` sample points from a subset of
points in `X`. Optimization is run for `gens` generations.
"""
function subLHCoptim(X,n::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)

    #preallocate memory
    nLarge, d = size(X)
    pop = [Array{Int}(undef,n,d) for i = 1:popsize+1]
    nextpop = similar(pop)
    fitness = Vector{Float64}(undef, popsize+1)
    fitnessInds = Vector{Int64}(undef, popsize+1)
    bestfits = Array{Float64}(undef, gens)
    dist = zeros(Float64,Int(n*(n-1)*0.5))

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
        subInds = sample(1:nLarge, n, replace = false)
        pop[i] = X[subInds,:]
        fitness[i] = AudzeEglaisObjective!(dist, pop[i])
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
            winnerInd = tournament!(fitnessInds,ntour,tour_inds,ptour)
            nextpop[i] = copy(pop[winnerInd])
        end
        

        #perform mutation
        for i = 2:popsize+1
            for j=1:muts[gens]
                subPoint = sample(1:n) #Randomly chosen point subLHC
                largePoint = sample(1:nLarge) #Randomly chosen point largeLHC

                while nextpop[i][subPoint,:] == X[largePoint,:]
                    largePoint = sample(1:nLarge) #Choose new random point in largeLHC
                end
                nextpop[i][subPoint,:] = X[largePoint,:]
            end
        end

        #evaluate fitness
        for i = 1:popsize+1
            fitness[i] = AudzeEglaisObjective!(dist, nextpop[i])
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
