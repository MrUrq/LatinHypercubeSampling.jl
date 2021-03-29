using LatinHypercubeSampling
using Random
using StableRNGs
using Test

@testset "AudzeEglais" begin

    LHC = [1 3; 3 1; 2 2]
    n = size(LHC,1)
    dist = zeros(Float64,Int(n*(n-1)*0.5))

    @test_logs (:warn,"AudzeEglaisObjective!(dist,LHC) is deprecated and does not differ from AudzeEglaisObjective(LHC)") AudzeEglaisObjective!(dist,LHC)
    @test AudzeEglaisObjective(LHC) ≈ 0.88888888888888888
end

@testset "mutateLHC" begin

    rng = StableRNGs.StableRNG(1)
    a = [0,0]
    LHC = randomLHC(rng,3,2)

    rng = StableRNGs.StableRNG(1)
    LatinHypercubeSampling.mutateLHC!(rng,LHC,a)
    @test LHC == [  2  1
                    3  2
                    1  3]
end

@testset "tournament!" begin

    pop = Array{Array{Float64}}(undef, 15)
    popfit = Array{Float64}(undef, 15)
    fitnessInds = Vector{Int64}(undef, 15)
    rng = StableRNGs.StableRNG(1)

    for i = 1:15
        LHC = randomLHC(rng,15,3)
        popfit[i] = AudzeEglaisObjective(LHC)
    end

    sortperm!(fitnessInds,popfit)

    tour_inds = Array{Int}(1:15)     #Storage of indices for tournament selection
    tour_fitinds = Array{Int}(1:15)  #Storage of fitness for tournament selection    
        
    @test LatinHypercubeSampling.tournament!(rng,fitnessInds,15,tour_inds,1) == 1
    @test LatinHypercubeSampling.tournament!(rng,fitnessInds,15,tour_inds,0) == 6
end

@testset "fixedcross" begin

    parone = [1,2,3,4,5,6,7,8]
    partwo = [4,2,5,1,6,8,3,7]
    offsprone = similar(parone)
    offsprtwo = similar(parone)
    randlock = ReentrantLock()
    rng = StableRNGs.StableRNG(1)
    @test [1,2,3,4,5,6,8,7] == LatinHypercubeSampling._fixedcross!(rng,offsprone,parone,partwo,randlock)
    rng = StableRNGs.StableRNG(1)
    @test [4,2,5,1,3,6,7,8]== LatinHypercubeSampling._fixedcross!(rng,offsprone,partwo,parone,randlock)

    rng = StableRNGs.StableRNG(1)
    @test [1,2,3,4,5,6,8,7] == LatinHypercubeSampling.fixedcross!(rng,offsprone,offsprtwo,parone,partwo,randlock)[1]
    rng = StableRNGs.StableRNG(1)
    @test [4,2,5,1,3,6,7,8] == LatinHypercubeSampling.fixedcross!(rng,offsprone,offsprtwo,parone,partwo,randlock)[2]
end

@testset "inversion" begin

    individual = [1,2,3,4]
    rng = StableRNGs.StableRNG(1)
    LatinHypercubeSampling.inversion!(rng,individual)

    @test [1,2,4,3] == individual

end

@testset "is optimization result still an LHC" for threading in (true, false)
    numPoints = 64
    numDims = 10
    numGens = 100

    X = LHCoptim(numPoints,numDims,numGens;popsize=20,ntour=2,ptour=0.8,periodic_ae=false,threading=threading)[1]
    n, d = size(X)
    for i = 1:d
        @test length(unique(X[:,i])) == n
    end

    X = LHCoptim(numPoints,numDims,numGens;popsize=20,ntour=2,ptour=0.8,periodic_ae=true,threading=threading)[1]
    n, d = size(X)
    for i = 1:d
        @test length(unique(X[:,i])) == n
    end

    X = subLHCoptim(X,numPoints÷2,numGens;popsize=20,ntour=2,ptour=0.8,periodic_ae=false)[1]
    n, d = size(X)
    for i = 1:d
        @test length(unique(X[:,i])) == n
    end

    X = subLHCoptim(X,numPoints÷2,numGens;popsize=20,ntour=2,ptour=0.8,periodic_ae=true)[1]
    n, d = size(X)
    for i = 1:d
        @test length(unique(X[:,i])) == n
    end

end

@testset "is the fittest individual kept" for threading in (true,false)
    numPoints = 64
    numDims = 10
    numGens = 100

    X,fitnesses = LHCoptim(numPoints,numDims,numGens;popsize=20,ntour=2,ptour=0.8,threading=threading)
    @test sort(fitnesses) == fitnesses
    fitnesses = subLHCoptim(X,numPoints÷2,numGens;popsize=20,ntour=2,ptour=0.8)[2]
    @test sort(fitnesses) == fitnesses
end

@testset "categorical LHC" for threading in (true, false)
    numPoints = 64
    numCat = 4
    W = 1.0
    catW= 0.01
    dims = [Continuous(),Continuous(),Categorical(numCat,catW)]
    
    numGens = 100

    X = LHCoptim!(randomLHC(numPoints,dims),numGens;dims=dims,interSampleWeight=W,threading=threading)[1]

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

@testset "scaleLHC" begin
    plan = randomLHC(5,2)
    scaled_plan = scaleLHC(plan,[(-1,1),(10,100)])    
    @test extrema(scaled_plan; dims=1) == [(-1.0, 1.0)  (10.0, 100.0)]
end

@testset "Global rng test" begin
    Random.seed!(1)
    before_lhc = rand(100)
    Random.seed!(1)
    
    LHCoptim(10,3,100,rng=Random.MersenneTwister(1))

    after_lhc = rand(100)

    @test before_lhc == after_lhc
end
