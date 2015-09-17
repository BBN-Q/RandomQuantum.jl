VERSION >= v"0.4.0-dev+6521" && __precompile__()
module RandomQuantum

using QuantumInfo

import Base.rand

export  # types
        GinibreEnsemble,
        FubiniStudyPureState,
        FubiniStudyMixedState,
        HilbertSchmidtMixedState,
        BuresMixedState,
        HaarUnitary,
        HaarCPTPMap,
        GUEUnitary,
        GUECPTPMap,
        # functions
        rand

"""
Type corresponding to the Ginibre distribution
of complex matrices.
"""
type GinibreEnsemble
    rows::Int64
    cols::Int64
end

GinibreEnsemble(dim::Number) = GinibreEnsemble(dim,dim)

function rand(dist::GinibreEnsemble)
    return (randn(dist.rows,dist.cols)+im*randn(dist.rows,dist.cols))/sqrt(2)
end

"""
Type corresponding to unitarily-invariant distribution
of pure states for some Hilbert space dimension.
"""
type FubiniStudyPureState
    dim::Int64
end

function rand(dist::FubiniStudyPureState)
    return normalize( rand(GinibreEnsemble(dist.size)) )
end

""" 
Type corresponding to distribution of mixed states obtained by
tracing out part of a pure state obtained from a FubiniStudy
distribution. When `dim == bath_dim`, this is identical to the
Hilbert-Schmidt distribution of mixed states. 
See, e.g., Zyczkowski et al., [J. Math. Phys. 52, 062201
(2011)](http://dx.doi.org/10.1063/1.3595693)
[arXiv:1010.3570](http://arxiv.org/abs/1010.3570).
"""
type FubiniStudyMixedState
    dim::Int64
    bath_dim::Int64
end

function rand(dist::FubiniStudyMixedState)
    return trace(projector(rand(FubiniStudyPureState(dist.dim+dist.bath_dim))),[dist.dim,dist.bath_dim],2)
end

"""
Type corresponding to distribution of mixed states according to the
Hilbert-Schmidt measure (induced by the Frobenius distance). This is
identical to the mixed state distribution induced by tracing out half
of a pure state distributed according to the Fubini-Study metric.
See, e.g., Zyczkowski et al., [J. Math. Phys. 52, 062201
(2011)](http://dx.doi.org/10.1063/1.3595693)
[arXiv:1010.3570](http://arxiv.org/abs/1010.3570).
"""
type HilbertSchmidtMixedState
    dim::Int64
end

function rand(dist::HilbertSchmidtMixedState)
    X = rand(GinibreEnsemble(dist.dim))
    M = X*X'
    return M/trace(M)
end

""" 
Type corresponding to distribution of mixed states according to the
Bures metric. 
See, e.g., Zyczkowski et al., [J. Math. Phys. 52, 062201
(2011)](http://dx.doi.org/10.1063/1.3595693)
[arXiv:1010.3570](http://arxiv.org/abs/1010.3570).
"""
type BuresMixedState
    dim::Int64
end

function rand(dist::BuresMixedState)
    G = rand(GinibreEnsemble(dist.dim))
    U = rand(HaarUnitary(dist.dim))
    H = (eye(dist.dim)+U)*G
    return H*H'/trace(H*H')
end

""" 
Type corresponding to the unitarily invariant distribution of unitary
transformations for some Hilbert space dimension.  See, e.g.,
Mezzadri, [Notices Amer Math Soc 54 4 592
(2007)](http://www.ams.org/notices/200705/fea-mezzadri-web.pdf).
"""
type HaarUnitary 
    dim::Int64 
end

function rand(dist::HaarUnitary)
    # TODO: replace with call to RandomMatrices?
    X = rand(GinibreEnsemble(dist.dim))
    Q,_ = qr(X,thin=false)
    d = diag(Q)
    d = d./abs(d)
    Q = Q./d
    return Q
end

"""
Type corresponding to the unitarily invariant distribution of
unitary transformations for system and bath, such that the bath
(initially in the ground state) is traced out.
See, e.g., Bruzda et al., [Phys. Lett. A 373, 320-324
(2009)](http://dx.doi.org/10.1016/j.physleta.2008.11.043),
[arXiv:0804.2361](http://arxiv.org/abs/0804.2361).
"""
type HaarCPTPMap
    dim::Int64
    bath_dim::Int64
    function HaarCPTPMap(dim,bath_dim)
        if bath_dim <= 1
            error("Bath must have dimension 2 or larger")
        end
    end
end

function rand(dist::HaarCPTPMap)
    X = rand(GinibreEnsemble(dist.dim^2,dist.bath_dim))
    W = X*X'
    Y = sqrtm(trace(W,[dist.dim,dist.dim],1))
    return (kron(eye(dist.dim),Y)\X)/kron(eye(dist.dim),Y)
end

"""
Type corresponding to the Gaussian Unitary Ensemble of complex
matrices.  
"""
type GUE 
    size::Int64
end

function rand(dist::GUE)
    X = rand(GinibreEnsemble(dist.size))
    return (X+X')/2
end

"""
Type corresponding to the integrated evolution of random Hamiltonians
(with unitarily invariant distribution) on some Hilbert space.
"""
type GUEUnitary
    dim::Int64
    α::Float64
end

# TODO: Instead of matrix exponentiation, could this be done by
#       sampling eigenvalues, exponetiating them, and then multiplying
#       by a Haar random unitary? If so, that should be cheaper
function rand(dist::GUEUnitary)
    rh = rand(GUE(dist.dim))
    return expm(-1im*dist.α*rh)
end

"""
Type corresponding to the distribution of completely positive trace
preserving maps obtained from the integrated evolution of random
Hamiltonians (with unitarily invariant distribution) on a
dimensional Hilbert space, such that bath (initially in the ground
state) is traced out after the evolution.
"""
type GUECPTPMap
    dim::Int64
    bath_dim::Int64
    α::Float64
end

function rand(GUECPTPMap)
    ru = rand(GUEUnitary(dist.dim+dist.bath_dim)) * kron(eye(d),ket(0,de))
    k = Matrix{Complex128}[]
    for ii=1:de
        push!(k, kron(eye(d),bra(ii-1,de)) * ru )
    end
    return kraus2liou(k)
end

end
