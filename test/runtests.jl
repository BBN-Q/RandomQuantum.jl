using QuantumInfo, RandomQuantum
using Base.Test

trnorm(A) = sum(svdvals(A))

#@testloop "Randomized GinUE tests" 
for _ in 1:10
    rows = rand(2:10)
    cols = rand(1:10)
    X = rand(GinibreEnsemble(rows,cols))
    @test var(vec(real(X))) > 0
    @test var(vec(imag(X))) > 0
end

#@testloop "Randomized FubiniStudyEnsemble tests" 
for _ in 1:10
    d = rand(2:10)
    ψ = rand(FubiniStudyPureState(d))
    @test var(real(ψ)) > 0
    @test var(imag(ψ)) > 0
    @test_approx_eq norm(ψ) 1.0
end

#@testloop "Randomized FSOE tests" 
for _ in 1:10
    d = rand(2:10)
    db = rand(2:10)
    ρ = rand(FubiniStudyMixedState(d,db))
    @test_approx_eq trnorm(ρ) 1
    @test ishermitian(ρ)
    @test ispossemidef(ρ)
    @test var(vec(real(ρ))) > 0
    @test var(vec(imag(ρ))) > 0
end

#@testloop "Randomized HSOE tests" 
for _ in 1:10
    d = rand(2:10)
    ρ = rand(HilbertSchmidtMixedState(d))
    @test_approx_eq trnorm(ρ) 1
    @test ishermitian(ρ)
    @test ispossemidef(ρ)
    @test var(vec(real(ρ))) > 0
    @test var(vec(imag(ρ))) > 0
end

#BuresMixedState
for _ in 1:10
    d = rand(2:10)
    ρ = rand(BuresMixedState(d))
    @test_approx_eq trnorm(ρ) 1
    @test ishermitian(ρ)
    @test ispossemidef(ρ)
    @test var(vec(real(ρ))) > 0
    @test var(vec(imag(ρ))) > 0
end

#ClosedHaarEnsemble
for _ in 1:10
    d = rand(2:10)
    U = rand(ClosedHaarEnsemble(d))
    @test_approx_eq_eps norm(U*U'-eye(d)) 0.0 1e-12
    @test_approx_eq_eps norm(U'*U-eye(d)) 0.0 1e-12
    @test_approx_eq_eps abs(eigvals(U)) ones(d) 1e-12
end

#OpenHaarEnsemble
for _ in 1:10
    d = 2 # rand(2:10)
    db = rand(2:10)
    E = rand(OpenHaarEnsemble(d,db))
    @test istp(E,tol=1e-14)
    @test iscp(E,tol=1e-14)
    @test ischannel(E,tol=1e-14)
end

#ClosedEvolution
for _ in 1:10
    α = rand()
    d = rand(2:10)
    U = rand(RandomClosedEvolution(d,α))
    @test_approx_eq_eps norm(U*U'-eye(d)) 0.0 1e-12
    @test_approx_eq_eps norm(U'*U-eye(d)) 0.0 1e-12
    @test_approx_eq_eps abs(eigvals(U)) ones(d) 1e-12
end

#OpenEvolution
for _ in 1:10
    α = rand()
    d = rand(2:10)
    db = rand(2:10)
    E = rand(RandomOpenEvolution(d,db,α))
    @test ischannel(E,tol=1e-12)
end

