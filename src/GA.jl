"""
    function mutateLHC!(LHC,mut_inds)
Perform one random mutation to the LHC. Interchange two randomly selected elements
from a randomly selected column.
"""
function mutateLHC!(LHC,mut_inds)

    n, d = size(LHC)         #Num points, num dimensions
    mut_range = 1:n          #Range of points per dim for mutation      
    mut_dim = 1:d            #Range of dimensions for randomly choosen dimension

    mut_inds =  sample!(mut_range, mut_inds, replace = false) #Two unique & random indices
    ind1,ind2 = mut_inds

    rnd_col = rand(mut_dim)   #Randomly choosen column (dimension)

    #Swap elements
    mut_tmp = LHC[ind1,rnd_col]
    LHC[ind1,rnd_col] = LHC[ind2,rnd_col]
    LHC[ind2,rnd_col] = mut_tmp

    return
end



"""
    function tournament!(pop,tour_num::Int,prob::Int)
Choose `tour_num` random, unique, individuals and choose the best individual
with probability `prob`. Return selected individual index
"""
function tournament!(popfit::Vector{Float64},tour_num,tour_inds,tour_fitinds,prob)

    tour_range = 1:length(popfit)

    #Unique individuals for tournament
    sample!(tour_range, tour_inds, replace = false)::Array{Int64,1}

    #Fitness among contenders sorted in tour_fitinds
    fitnesses = view(popfit,tour_inds)
    sortperm!(tour_fitinds,fitnesses,rev=true)

    #Return the selected index of the individual
    for i = 1:tour_num-1
        if rand() <= prob
            fittest = tour_fitinds[i]
            return tour_inds[fittest]    #Winner!
        end
    end
    return tour_inds[tour_fitinds[end]]  #Lucky chum!

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
function _fixedcross(parone,partwo)

    #initialise offspring
    offspr = zeros(typeof(parone[1]),length(parone))

    #generate a random location in the gene
    n = length(parone)
    loc = sample(1:n-1)

    offspr[1:loc] = parone[1:loc]
    i = loc+1
    while i < n+1
        for j = 1:n
            @compat x = coalesce(findfirst(isequal(partwo[j]), offspr), 0)
            if (x == 0) && (offspr[i] == 0)
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
function fixedcross(parone,partwo)

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
    tmp = reverse(individual[invRange], 1)    
    individual[invRange] = tmp

    return

end
