module RandomQuantum

using  QuantumInfo
import QuantumInfo: eye
using LinearAlgebra: tr, eigen, normalize, qr, diag, diagm
import Base.rand

export  # types
        GinibreEnsemble,
        FubiniStudyPureState,
        FubiniStudyMixedState,
        HilbertSchmidtMixedState,
        BuresMixedState,
        ClosedHaarEnsemble,
        OpenHaarEnsemble,
        RandomClosedEvolution,
        RandomOpenEvolution,
        # functions
        rand

function eig(m)
    r = eigen(m)
    return (r.values, r.vectors)
end

"""
RandomQuantum.GinibreEnsemble(rows,cols)

Type corresponding to the Ginibre distribution
of complex matrices.
"""
struct GinibreEnsemble
    rows::Int64
    cols::Int64
end

GinibreEnsemble(dim::Number) = GinibreEnsemble(dim,dim)

function rand(dist::GinibreEnsemble)
    return (randn(dist.rows,dist.cols)+im*randn(dist.rows,dist.cols))/sqrt(2)
end

"""
RandomQuantum.FubiniStudyPureState(dim)

Type corresponding to unitarily-invariant distribution
of pure states for some Hilbert space dimension.
"""
struct FubiniStudyPureState
    dim::Int64
end

function rand(dist::FubiniStudyPureState)
    return normalize( vec(rand(GinibreEnsemble(dist.dim,1))) )
end

"""
RandomQuantum.FubiniStudyMixedState(dim, bath_dim)

Type corresponding to distribution of mixed states obtained by
tracing out part of a pure state obtained from a FubiniStudy
distribution. When `dim == bath_dim`, this is identical to the
Hilbert-Schmidt distribution of mixed states.
See, e.g., Zyczkowski et al., [J. Math. Phys. 52, 062201
(2011)](http://dx.doi.org/10.1063/1.3595693)
[arXiv:1010.3570](http://arxiv.org/abs/1010.3570).
"""
struct FubiniStudyMixedState
    dim::Int64
    bath_dim::Int64
end

function rand(dist::FubiniStudyMixedState)
    ψ = rand(FubiniStudyPureState(dist.dim*dist.bath_dim))
    return partialtrace(projector(ψ),[dist.dim,dist.bath_dim],2)
end

"""
RandomQuantum.HilbertSchmidtState(dim)

Type corresponding to distribution of mixed states according to the
Hilbert-Schmidt measure (induced by the Frobenius distance). This is
identical to the mixed state distribution induced by tracing out half
of a pure state distributed according to the Fubini-Study metric.
See, e.g., Zyczkowski et al., [J. Math. Phys. 52, 062201
(2011)](http://dx.doi.org/10.1063/1.3595693)
[arXiv:1010.3570](http://arxiv.org/abs/1010.3570).
"""
struct HilbertSchmidtMixedState
    dim::Int64
end

function rand(dist::HilbertSchmidtMixedState)
    X = rand(GinibreEnsemble(dist.dim))
    M = X*X'
    return M/tr(M)
end

"""
RandomQuantum.BuresMixedState(dim)

Type corresponding to distribution of mixed states according to the
Bures metric.
See, e.g., Zyczkowski et al., [J. Math. Phys. 52, 062201
(2011)](http://dx.doi.org/10.1063/1.3595693)
[arXiv:1010.3570](http://arxiv.org/abs/1010.3570).
"""
struct BuresMixedState
    dim::Int64
end

function rand(dist::BuresMixedState)
    G = rand(GinibreEnsemble(dist.dim))
    U = rand(ClosedHaarEnsemble(dist.dim))
    H = (eye(dist.dim)+U)*G
    return H*H'/tr(H*H')
end

"""
RandomQuantum.ClosedHaarEnsemble(dim)

Type corresponding to the unitarily invariant distribution of unitary
transformations for some Hilbert space dimension.  See, e.g.,
Mezzadri, [Notices Amer Math Soc 54 4 592
(2007)](http://www.ams.org/notices/200705/fea-mezzadri-web.pdf).
"""
struct ClosedHaarEnsemble
    dim::Int64
end

function rand(dist::ClosedHaarEnsemble)
    X = rand(GinibreEnsemble(dist.dim))
    Q,_ = qr(X)
    d = diag(Q)
    d = d ./ abs.(d)
    Q = Q ./ d
    return Q
end

"""
RandomQuantum.OpenHaarEnsemble(dim, bath_dim)

Type corresponding to the unitarily invariant distribution of
unitary transformations for system and bath, such that the bath
(initially in the ground state) is traced out.
See, e.g., Bruzda et al., [Phys. Lett. A 373, 320-324
(2009)](http://dx.doi.org/10.1016/j.physleta.2008.11.043),
[arXiv:0804.2361](http://arxiv.org/abs/0804.2361).
"""
struct OpenHaarEnsemble
    dim::Int64
    bath_dim::Int64
end

function rand(dist::OpenHaarEnsemble)
    X = rand(GinibreEnsemble(dist.dim^2,dist.bath_dim))
    W = X*X'
    W = W/tr(W)
    #Y = sqrt(tr(W,[dist.dim,dist.dim],1))
    #IY = kron(eye(dist.dim),Y)
    Y = sqrt(partialtrace(W,[dist.dim,dist.dim],2))
    IY = kron(Y,eye(dist.dim))
    R = IY\W/IY
    return choi2liou(R)/dist.dim
end

"""
RandomQuantum.GUE(dim)

Type corresponding to the Gaussian Unitary Ensemble of complex
matrices.
"""
struct GUE
    size::Int64
end

function rand(dist::GUE)
    X = rand(GinibreEnsemble(dist.size))
    return (X+X')/2
end

"""
RandomQuantum.RandomClosedEvolution(dim, α)

Type corresponding to the integrated evolution of random Hamiltonians
(with unitarily invariant distribution) on some Hilbert space.
"""
struct RandomClosedEvolution
    dim::Int64
    α::Float64
end

# TODO: Instead of matrix exponentiation, could this be done by
#       sampling eigenvalues, exponetiating them, and then multiplying
#       by a Haar random unitary? If so, that should be cheaper
function rand(dist::RandomClosedEvolution)
    evals,U = eig(rand(GUE(dist.dim)))
    return U*diagm(0 => exp.(-1im*dist.α*evals))*U'
end

"""
RandomQuantum.RandomOpenEvolution(dim, bath_dim, α)

Type corresponding to the distribution of completely positive trace
preserving maps obtained from the integrated evolution of random
Hamiltonians (with unitarily invariant distribution) on a
dimensional Hilbert space, such that bath (initially in the ground
state) is traced out after the evolution.
"""
struct RandomOpenEvolution
    dim::Int64
    bath_dim::Int64
    α::Float64
end

function rand(dist::RandomOpenEvolution)
    d = dist.dim
    de = dist.bath_dim
    ru = rand(RandomClosedEvolution(d*de,dist.α)) * kron(eye(d),ket(0,de))
    k = Matrix{ComplexF64}[]
    for ii=1:de
        push!(k, kron(eye(d), bra(ii-1,de)) * ru )
    end
    return kraus2liou(k)
end

end
