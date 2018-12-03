using QuantumInfo, RandomQuantum
using LinearAlgebra: svdvals, norm, ishermitian, eigvals
using StatsBase: var

using Test

trnorm(A) = sum(svdvals(A))

@testset "Randomized GinUE tests" begin
for _ in 1:10
    rows = rand(2:10)
    cols = rand(1:10)
    X = rand(GinibreEnsemble(rows,cols))
    @test var(vec(real(X))) > 0
    @test var(vec(imag(X))) > 0
end
end

@testset "Randomized FubiniStudyEnsemble tests" begin
for _ in 1:10
    d = rand(2:10)
    ψ = rand(FubiniStudyPureState(d))
    @test var(real(ψ)) > 0
    @test var(imag(ψ)) > 0
    @test isapprox(norm(ψ), 1.0)
end
end

@testset "Randomized FSOE tests" begin
for _ in 1:10
    d = rand(2:10)
    db = rand(2:10)
    ρ = rand(FubiniStudyMixedState(d,db))
    @test isapprox(trnorm(ρ), 1)
    @test ishermitian(ρ)
    @test ispossemidef(ρ)
    @test var(vec(real(ρ))) > 0
    @test var(vec(imag(ρ))) > 0
end
end

@testset "Randomized HSOE tests" begin
for _ in 1:10
    d = rand(2:10)
    ρ = rand(HilbertSchmidtMixedState(d))
    @test isapprox(trnorm(ρ), 1)
    @test ishermitian(ρ)
    @test ispossemidef(ρ)
    @test var(vec(real(ρ))) > 0
    @test var(vec(imag(ρ))) > 0
end
end

@testset "BuresMixedState" begin
for _ in 1:10
    d = rand(2:10)
    ρ = rand(BuresMixedState(d))
    @test isapprox(trnorm(ρ), 1)
    @test ishermitian(ρ)
    @test ispossemidef(ρ)
    @test var(vec(real(ρ))) > 0
    @test var(vec(imag(ρ))) > 0
end
end

@testset "ClosedHaarEnsemble" begin
for _ in 1:10
    d = rand(2:10)
    U = rand(ClosedHaarEnsemble(d))
    @test isapprox(U*U', QuantumInfo.eye(d), atol=1e-12)
    @test isapprox(U'*U, QuantumInfo.eye(d), atol=1e-12)
    @test isapprox(abs.(eigvals(U)),  ones(d), atol=1e-12)
end
end

@testset "OpenHaarEnsemble" begin
for _ in 1:10
    d = 2 # rand(2:10)
    db = rand(2:10)
    E = rand(OpenHaarEnsemble(d,db))
    @test istp(E,tol=1e-14)
    @test iscp(E,tol=1e-14)
    @test ischannel(E,tol=1e-14)
end
end

@testset "ClosedEvolution" begin
for _ in 1:10
    α = rand()
    d = rand(2:10)
    U = rand(RandomClosedEvolution(d,α))
    @test isapprox(U*U', QuantumInfo.eye(d), atol=1e-12)
    @test isapprox(U'*U, QuantumInfo.eye(d), atol=1e-12)
    @test isapprox(abs.(eigvals(U)), ones(d), atol=1e-12)
end
end

@testset "OpenEvolution" begin
for _ in 1:10
    α = rand()
    d = rand(2:10)
    db = rand(2:10)
    E = rand(RandomOpenEvolution(d,db,α))
    @test ischannel(E,tol=1e-12)
end
end
