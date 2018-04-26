using LatinHypercubeSampling
using Base.Test

@testset "randomLHC" begin

    srand(1)
    @test LHC = randomLHC(3,2) == [1  3
                                   3  1
                                   2  2]


end

@testset "AudzeEglias" begin

    srand(1)
    LHC = randomLHC(3,2)
    n = size(LHC,1)
    dist = zeros(Float64,Int(n*(n-1)*0.5))

    @test AudzeEgliasObjective(dist,LHC) ≈ 0.88888888888888888

end

@testset "mutateLHC" begin

    srand(1)
    LHC = randomLHC(3,2)
    mutateLHC!(LHC)
    @test LHC == [3  3
                  1  1
                  2  2]
end

@testset "tournament" begin

    pop = Array{Array{Float64}}(15)
    popfit = Array{Float64}(15)
    srand(1)

    for i = 1:15
        LHC = randomLHC(15,3)
        n = size(LHC,1)
        dist = zeros(Float64,Int(n*(n-1)*0.5))
        popfit[i] = AudzeEgliasObjective(dist,LHC)
    end

    @test tournament(popfit,15,1) == 11
    @test tournament(popfit,15,0) == 7
end

@testset "cyclecross" begin

    parone = [3,4,8,2,7,1,6,5]
    partwo = [4,2,5,1,6,8,3,7]

    @test parone == _cyclecross(parone,partwo)
    @test partwo == _cyclecross(partwo,parone)
    @test parone == cyclecross(parone,partwo)[1]
    @test partwo == cyclecross(parone,partwo)[2]
end

@testset "fixedcross" begin

    parone = [1,2,3,4,5,6,7,8]
    partwo = [4,2,5,1,6,8,3,7]
    srand(2)
    @test [1,2,4,5,6,8,3,7] == _fixedcross(parone,partwo)
    srand(2)
    @test [4,2,1,3,5,6,7,8] == _fixedcross(partwo,parone)

    srand(2)
    @test [1,2,4,5,6,8,3,7] == fixedcross(parone,partwo)[1]
    srand(2)
    @test [4,2,1,3,5,6,7,8] == fixedcross(parone,partwo)[2]
end


@testset "inversion" begin

    individual = [1,2,3,4,5,6,7,8]
    srand(1)
    inversion!(individual)

    @test [1,2,3,7,6,5,4,8] == individual

end


@testset "is optimization result still an LHC" begin
    numPoints = 64
    numDims = 10
    numGens = 100

    X = LHCoptim(numPoints,numDims,numGens;popsize=20,ntour=2,ptour=0.8)[1]
    n, d = size(X)
    for i = 1:d
        @test length(unique(X[:,1])) == n
    end

    X = subLHCoptim(X,numDims,numGens;popsize=20,ntour=2,ptour=0.8)[1]
    n, d = size(X)
    for i = 1:d
        @test length(unique(X[:,1])) == n
    end

    X = refineLHCoptim(X,numGens;popsize=20)[1]
    n, d = size(X)
    for i = 1:d
        @test length(unique(X[:,1])) == n
    end
end
