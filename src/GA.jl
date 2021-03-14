"""
    function mutateLHC!(rng,LHC,mut_inds)
Perform one random mutation to the LHC. Interchange two randomly selected elements
from a randomly selected column.
"""
function mutateLHC!(rng,LHC,mut_inds)

    n, d = size(LHC)         #Num points, num dimensions
    mut_range = 1:n          #Range of points per dim for mutation      
    mut_dim = 1:d            #Range of dimensions for randomly choosen dimension

    mut_inds = sample!(rng, mut_range, mut_inds, replace = false) #Two unique & random indices
    ind1,ind2 = mut_inds

    rnd_col = rand(rng,mut_dim)   #Randomly choosen column (dimension)

    #Swap elements
    mut_tmp = LHC[ind1,rnd_col]
    LHC[ind1,rnd_col] = LHC[ind2,rnd_col]
    LHC[ind2,rnd_col] = mut_tmp

    return
end

"""
    function tournament!(rng,popfit_inds,tour_num,tour_inds,prob)
Choose `tour_num` random, unique, individuals and choose the best individual
with probability `prob`. Return selected individual index
"""
function tournament!(rng,popfit_inds,tour_num,tour_inds,prob)

    sample!(rng, popfit_inds, tour_inds, replace = false)::Array{Int64,1}
    
    sort!(tour_inds, rev=true)

    for i in 1:tour_num-1
        if rand(rng) <= prob
            return popfit_inds[tour_inds[i]]     #Winner!
        end
    end
    return popfit_inds[tour_inds[end]]           #Lucky chum!
end

"""
    function _fixedcross!(rng,offspr,parone,partwo)
Fixed point crossover of two parents to create one offspring.
"""
function _fixedcross!(rng,offspr,parone,partwo,randlock)

    offspr .= 0

    #generate a random location in the gene
    n = length(parone)
    lock(randlock)
        loc = sample(rng,1:n-1)
    unlock(randlock)
    offspr[1:loc] = parone[1:loc]
    i = loc+1
    while i < n+1
        for j = 1:n
            x = findfirst(isequal(partwo[j]), offspr)
            if (x isa Nothing) && (offspr[i] == 0)
                offspr[i] = partwo[j]
                i += 1
            end
        end
    end

    return offspr

end

"""
    function fixedcross!(rng,offsprone,offsprtwo,parone,partwo)
Fixed point crossover of two parents to create two offspring.
"""
function fixedcross!(rng,offsprone,offsprtwo,parone,partwo,randlock=ReentrantLock())

    _fixedcross!(rng,offsprone,parone,partwo,randlock)
    _fixedcross!(rng,offsprtwo,partwo,parone,randlock)

    return offsprone, offsprtwo

end

"""
    function inversion!(rng,individual)
Reverse the values in a gene between two random cut-off points.
"""
function inversion!(rng,individual)

    #generate two, unique, random locations in the gene
    loc = sample(rng, 1:length(individual), 2, replace = false)

    #flip the values in the range in place
    reverse!(individual, minimum(loc),maximum(loc))    

    return

end
