module LatinHypercubeSampling

export  randomLHC,
        AudzeEgliasObjective,
        LHCoptim,
        LHCoptim!,
        subLHCoptim,
        subLHCindex,
        refineLHCoptim,



include("GA.jl")


using StatsBase

"""
    function randomLHC(n::Int, d::Int)
Generate a random Latin Hypercube with `d` dimensions and `n` sample points.
"""
function randomLHC(n::Int,d::Int)

    LHC = Array{Int}(n,d)
    for i = 1:d
        LHC[:,i] = randperm(n)
    end

    return LHC

end



"""
    function AudzeEgliasObjective(dist,LHC::Array)
Return the scalar which should be maximized when using the Audze-Eglias
distance as the objective function. Note this is the inverse of the typical
Audze-Eglias distance which normally is minimized.
"""
function AudzeEgliasObjective(dist,LHC::Array)

    n = size(LHC,1) #Num points
    d = size(LHC,2) #Num dimensions
    dist .= 0.0

    # Squared l-2 norm of distances between all (unique) points
    l = 0
    for i = 2:n
        for j = 1:i-1
            l += 1
            for k = 1:d
                dist[l] += (LHC[i,k]-LHC[j,k])^2
            end
            dist[l] = 1/dist[l]
        end
    end

    output = 1/(sum(dist))

end



"""
    function AudzeEgliasObjective(LHC::Array)
Same as AudzeEgliasObjective(dist,LHC::Array) but creating a new distance array.
"""
function AudzeEgliasObjective(LHC::Array)

    n = size(LHC,1) #Num points
    dist = zeros(Float64,Int(n*(n-1)*0.5))
    AudzeEgliasObjective(dist,LHC::Array)

end



"""
    function LHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)
Produce an optimized Latin Hyper Cube with `d` dimensions and `n` sample points.
Optimization is run for `gens` generations.
"""
function LHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)

    #populate first individual
    X = randomLHC(n,d)

    LHCoptim!(X,gens,popsize=popsize,ntour=ntour,ptour=ptour)

end



"""
    function LHCoptim!(X::Array{Int,2},n::Int,d::Int,gens;ntour::Int=2,ptour=0.8)
Same as LHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8) but using an
existing population. Useful for continued optimization.
"""
function LHCoptim!(X::Array{Int,2},gens;popsize=100,ntour::Int=2,ptour=0.8)

    #preallocate memory
    n,d = size(X)
    pop = Array{Int}(popsize,n,d)

    pop[1,:,:] = X
    for i = 2:popsize
        pop[i,:,:] = randomLHC(n,d)
    end


    #preallocate memory
    nextpop = Array{Int}(popsize,n,d)
    fitness = Vector{Float64}(popsize)
    bestfits = Array{Float64}(gens)
    dist = zeros(Float64,Int(n*(n-1)*0.5))


    #crossover for even population size
    popEven = iseven(popsize)*-1


    #dynamic mutation rate table
    muts = Array{Int}(gens)
    for l = 1:gens
        muts[l] = round(Int,-n/(gens*0.75)*l+n)
        if muts[l] < 1
            muts[l] = 1
        end
    end


    #populate first populations evaluate fitness
    for i = 1:popsize
        fitness[i] = AudzeEgliasObjective(dist, pop[i,:,:])
    end


    #save the best individual and it's fitness
    bestfit, bestind = findmax(fitness)
    nextpop[1,:,:] = pop[bestind,:,:]
    bestfits[1] = bestfit


    #iterate for gens generations
    for k = 1:gens

        #tournament selection
        for i = 2:popsize
            winnerInd = tournament(fitness,ntour,ptour)
            nextpop[i,:,:] = pop[winnerInd,:,:]
        end

        #create children from crossover
        for i = 2:2:popsize+popEven
            for j = 1:d
                if rand() < 1/d
                    nextpop[i,:,j], nextpop[i+1,:,j] = fixedcross(pop[i,:,j],pop[i+1,:,j])
                end
            end
        end

        #perform inversion and mutation
        for i = 2:popsize
            for j = 1:d
                if rand() < 1/d     #probability for inversion
                    inversion!(view(nextpop,i,:,j))
                end
            end

            for j=1:muts[gens]
                mutateLHC!(view(nextpop,i,:,:))
            end

        end

        #evaluate fitness
        for i = 1:popsize
            fitness[i] = AudzeEgliasObjective(dist, nextpop[i,:,:])
        end

        #set altered population to current
        pop = copy(nextpop)

        #set the first individual to the best and save the fitness
        bestfit, bestind = findmax(fitness)
        nextpop[1,:,:] = pop[bestind,:,:]
        bestfits[k] = bestfit
    end

    return nextpop[1,:,:], bestfits

end



"""
    function subLHCoptim(X,n::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)
Produce an optimized Latin Hyper Cube with `n` sample points from a subset of
points in `X`. Optimization is run for `gens` generations.
"""
function subLHCoptim(X,n::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)

    #preallocate memory
    nLarge, d = size(X)
    pop = Array{Int}(popsize+1,n,d)
    nextpop = Array{Int}(popsize+1,n,d)
    fitness = Vector{Float64}(popsize+1)
    bestfits = Array{Float64}(gens)
    dist = zeros(Float64,Int(n*(n-1)*0.5))


    #dynamic mutation rate table
    muts = Array{Int}(gens)
    for l = 1:gens
        muts[l] = round(Int,-n/(gens*0.75)*l+n)
        if muts[l] < 1
            muts[l] = 1
        end
    end


    #populate first population with random subLHC's and evaluate fitness
    for i = 1:popsize+1
        subInds = sample(1:nLarge, n, replace = false)
        pop[i,:,:] = X[subInds,:]
        fitness[i] = AudzeEgliasObjective(dist, pop[i,:,:])
    end


    #save the best individual and it's fitness
    bestfit, bestind = findmax(fitness)
    nextpop[1,:,:] = pop[bestind,:,:]
    bestfits[1] = bestfit


    #iterate for gens generations
    for k = 1:gens

        #tournament selection
        for i = 2:popsize+1
            winnerInd = tournament(fitness,ntour,ptour)
            nextpop[i,:,:] = pop[winnerInd,:,:]
        end

        #perform mutation
        for i = 2:popsize+1
            for j=1:muts[gens]
                subPoint = sample(1:n) #Randomly chosen point subLHC
                largePoint = sample(1:nLarge) #Randomly chosen point largeLHC

                while nextpop[i,subPoint,:] == X[largePoint,:]
                    largePoint = sample(1:nLarge) #Choose new random point in largeLHC
                end
                nextpop[i,subPoint,:] = X[largePoint,:]
            end
        end

        #evaluate fitness
        for i = 1:popsize+1
            fitness[i] = AudzeEgliasObjective(dist, nextpop[i,:,:])
        end

        #set altered population to current
        pop = copy(nextpop)

        #set the first individual to the best and save the fitness
        bestfit, bestind = findmax(fitness)
        nextpop[1,:,:] = pop[bestind,:,:]
        bestfits[k] = bestfit
    end

    return nextpop[1,:,:], bestfits

end



"""
    function subLHCindex(X,Xsub)
Index in the large LHC to get the subLHC.
"""
function subLHCindex(X,Xsub)

    nsub = size(Xsub,1)
    subInds = Array{Int}(nsub)

    for i = 1:nsub
        subInds[i] = find(all(Xsub[i,:]' .== X, 2))[1]
    end

    return subInds
end



"""
    function refineLHCoptim(X,gens;popsize::Int=100)
Refine an existing plan by mutation only.
"""
function refineLHCoptim(X,gens;popsize::Int=100)

    n, d = size(X)

    #preallocate memory
    pop = Array{Int}(popsize,n,d)
    nextpop = Array{Int}(popsize,n,d)
    fitness = Array{Float64}(popsize)
    bestfits = Array{Float64}(gens)
    dist = zeros(Float64,Int(n*(n-1)*0.5))

    initFit = AudzeEgliasObjective(dist, X)
    for i = 1:popsize
        nextpop[i,:,:] = X
        fitness[i] = initFit
    end


    #dynamic mutation rate table
    muts = Array{Int}(gens)
    for l = 1:gens
        muts[l] = round(Int,-n/(gens*0.75)*l+n)
        if muts[l] < 1
            muts[l] = 1
        end
    end


    #iterate for gens generations
    for k = 1:gens

        #perform mutation and evaluate fitness
        for i = 2:popsize
            for j=1:muts[gens]
                mutateLHC!(view(nextpop,i,:,:))
            end
        end
        for i = 1:popsize
            fitness[i] = AudzeEgliasObjective(dist, nextpop[i,:,:])
        end

        #set altered population to current
        pop = copy(nextpop)

        #set the population to the best individual if the individuals fitness is
        #below the threshold of the best
        bestfit, bestind = findmax(fitness)
        nextpop[1,:,:] = pop[bestind,:,:]
        for i = 1:popsize
            if fitness[i] < 0.98*bestfit
                nextpop[i,:,:] = pop[bestind,:,:]
            end
        end
        bestfits[k] = bestfit
    end

    return nextpop[1,:,:], bestfits

end


end # module
