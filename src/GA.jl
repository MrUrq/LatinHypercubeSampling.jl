"""
    function mutateLHC!(LHC)
Perform one random mutation to the LHC. Interchange two randomly selected elements
from a randomly selected column.
"""
function mutateLHC!(LHC)

    n, d = size(LHC) #Num points, num dimensions

    col = sample(1:d) #Randomly choosen column
    ind1 ,ind2 =  sample(1:n, 2, replace = false) #Two unique & random indices

    #Swap elements
    tmp = LHC[ind1,col]
    LHC[ind1,col] = LHC[ind2,col]
    LHC[ind2,col] = tmp

    return

end



"""
    function tournament(pop,numtour::Int,prob::Int)
Choose `numtour` random, unique, individuals and choose the best individual
with probability `prob`. Return selected individual index
"""
function tournament(popfit::Vector{Float64},numtour::Int,prob)

    p = length(popfit)

    #Unique individuals for tournament
    tourInds = sample(1:p, numtour, replace = false)

    #Fitness among contenders sorted in fitInds
    fitnesses = popfit[tourInds]
    fitInds = sortperm(fitnesses,rev=true)

    #Return the selected index of the individual
    for i = 1:numtour-1
        if rand() <= prob
            fittest = fitInds[i]
            return tourInds[fittest]    #Winner!
        end
    end
    return tourInds[fitInds[end]]       #Lucky chum!

end



"""
    function _cyclecross(parone,partwo)
Cyclecrossover of two parents to create one offspring.
"""
function _cyclecross(parone::Vector,partwo::Vector)

    #initialise offspring
    @compat offspr = Vector{typeof(parone[1])}(undef, length(parone))
    @compat visited = BitArray(undef, length(parone)).=false

    #first value is direct copy of the first parent
    offspr[1] = parone[1]
    ind = 1

    #loop over values until all possible visits are made
    while visited[ind] == false
      visited[ind] = true
      @compat ind = coalesce(findfirst(isequal(partwo[ind]), parone), 0)
      offspr[ind] = parone[ind]
    end

    #use remaining values from the second parent
    flipbits!(visited)
    offspr[visited] .= partwo[visited]

    return offspr

end



"""
    function cyclecross(parone,partwo)
Cyclecrossover of two parents to create two offspring.
"""
function cyclecross(parone::Vector,partwo::Vector)

    offsprone = _cyclecross(parone,partwo)
    offsprtwo = _cyclecross(partwo,parone)

    return offsprone, offsprtwo

end



"""
    function _fixedcross(parone,partwo)
Fixed point crossover of two parents to create one offspring.
"""
function _fixedcross(parone::Vector,partwo::Vector)

    #initialise offspring
    offspr = zeros(typeof(parone[1]),length(parone))

    #generate a random location in the gene
    n = length(parone)
    loc = sample(1:n-1)

    offspr[1:loc] = parone[1:loc]
    i = loc+1
    while i < n+1
        for j = 1:n
            @compat x =  coalesce(findfirst(isequal(partwo[j]), offspr), 0)            
            if x == 0 && offspr[i] == 0
                offspr[i] = partwo[j]
                i += 1
            end
        end
    end

    return offspr

end



"""
    function fixedcross(parone,partwo)
Fixed point crossover of two parents to create two offspring.
"""
function fixedcross(parone::Vector,partwo::Vector)

    offsprone = _fixedcross(parone,partwo)
    offsprtwo = _fixedcross(partwo,parone)

    return offsprone, offsprtwo

end



"""
    function inversion!(individual)
Reverse the values in a gene between two random cut-off points.
"""
function inversion!(individual)

    #generate two, unique, random locations in the gene
    n = length(individual)
    invLocs = sample(1:n, 2, replace = false)
    sort!(invLocs)

    #range of inversion
    invRange = invLocs[1]:invLocs[2]

    #flip the values in the range and assign them to the vector in place
    tmp = flipdim(individual[invRange],1)
    individual[invRange] = tmp

    return

end
