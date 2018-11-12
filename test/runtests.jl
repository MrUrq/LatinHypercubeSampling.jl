using LatinHypercubeSampling
using Random
using Test


@testset "randomLHC" begin

    Random.seed!(1)
    @test LHC = randomLHC(3,2) == [1  3
                                   3  1
                                   2  2]


end

@testset "AudzeEglais" begin

    Random.seed!(1)
    LHC = randomLHC(3,2)
    n = size(LHC,1)
    dist = zeros(Float64,Int(n*(n-1)*0.5))

    @test AudzeEglaisObjective!(dist,LHC) ≈ 0.88888888888888888
    @test AudzeEglaisObjective(LHC) ≈ 0.88888888888888888
end

@testset "mutateLHC" begin

    a = [0,0]
    Random.seed!(1)
    LHC = randomLHC(3,2)
    Random.seed!(1)
    LatinHypercubeSampling.mutateLHC!(LHC,a)
    @test LHC == [  1  2
                    3  1
                    2  3]
end

@testset "tournament!" begin

    pop = Array{Array{Float64}}(undef, 15)
    popfit = Array{Float64}(undef, 15)
    fitnessInds = Vector{Int64}(undef, 15)
    Random.seed!(1)

    for i = 1:15
        LHC = randomLHC(15,3)
        n = size(LHC,1)
        dist = zeros(Float64,Int(n*(n-1)*0.5))
        popfit[i] = AudzeEglaisObjective!(dist,LHC)
    end

    sortperm!(fitnessInds,popfit)

    tour_inds = Array{Int}(1:15)     #Storage of indices for tournament selection
    tour_fitinds = Array{Int}(1:15)  #Storage of fitness for tournament selection    
        
    @test LatinHypercubeSampling.tournament!(fitnessInds,15,tour_inds,1) == 11
    @test LatinHypercubeSampling.tournament!(fitnessInds,15,tour_inds,0) == 7
end


@testset "fixedcross" begin

    parone = [1,2,3,4,5,6,7,8]
    partwo = [4,2,5,1,6,8,3,7]
    offsprone = similar(parone)
    offsprtwo = similar(parone)
    Random.seed!(1)
    @test [1,2,3,4,5,6,7,8] == LatinHypercubeSampling._fixedcross!(offsprone,parone,partwo)
    Random.seed!(1)
    @test [4,2,5,1,6,8,3,7] == LatinHypercubeSampling._fixedcross!(offsprone,partwo,parone)

    Random.seed!(1)
    @test [1,2,3,4,5,6,7,8] == LatinHypercubeSampling.fixedcross!(offsprone,offsprtwo,parone,partwo)[1]
    Random.seed!(1)
    @test [4,2,5,1,3,6,7,8] == LatinHypercubeSampling.fixedcross!(offsprone,offsprtwo,parone,partwo)[2]
end


@testset "inversion" begin

    individual = [1,2,3,4]
    Random.seed!(1)
    LatinHypercubeSampling.inversion!(individual)

    @test [1,2,4,3] == individual

end


@testset "is optimization result still an LHC" begin
    numPoints = 64
    numDims = 10
    numGens = 100

    X = LHCoptim(numPoints,numDims,numGens;popsize=20,ntour=2,ptour=0.8)[1]
    n, d = size(X)
    for i = 1:d
        @test length(unique(X[:,i])) == n
    end

    X = subLHCoptim(X,numPoints÷2,numGens;popsize=20,ntour=2,ptour=0.8)[1]
    n, d = size(X)
    for i = 1:d
        @test length(unique(X[:,i])) == n
    end

end

@testset "is the fittest individual kept" begin
    numPoints = 64
    numDims = 10
    numGens = 100

    X,fitnesses = LHCoptim(numPoints,numDims,numGens;popsize=20,ntour=2,ptour=0.8)
    @test sort(fitnesses) == fitnesses
    fitnesses = subLHCoptim(X,numPoints÷2,numGens;popsize=20,ntour=2,ptour=0.8)[2]
    @test sort(fitnesses) == fitnesses
end


@testset "categorical LHC" begin
    numPoints = 64
    numCat = 4
    W = 1.0
    catW= 0.01
    dims = [Continuous(),Continuous(),Categorical(numCat,catW)]
    
    numGens = 100

    X = LHCoptim!(randomLHC(numPoints,dims),numGens;dims=dims,interSampleWeight=W)[1]

    n, d = size(X)
    for (i, dim) in enumerate(dims)
        if typeof(dim) == Continuous
            @test length(unique(X[:,i])) == n
        else
            @test length(unique(X[:,i])) == numCat
        end
    end

    @test Categorical(1) == Categorical(1, 0.0)

end



@testset "subLHCindex" begin
    X = randomLHC(5,2)
    Xsub = X[1:2:end,:]
    subInds = subLHCindex(X,Xsub)
    @test subInds == [1,3,5]
end
