module RandomQuantum

using QuantumInfo,
      RandomMatrices

import Base.rand

export # types
       FubuniStudyState,
       FubiniStudyMixedState,
       HilbertSchmidtState,
       HaarUnitary,
       HaarCPTPMap,
       GUEUnitary,
       GUECPTPMap,
       # functions
       rand

"""
Type corresponding to unitarily-invariant distribution
of pure states for some Hilbert space dimension.
"""
type FubiniStudyPureState
    dim::Int64
end

function rand(dist::FubiniStudyPureState)
  return normalize( randn(dist.dim)+ 1im*randn(dist.dim) )
end

function rand(dist::FubiniStudyPureState, n::Int)
  if n<1
      error("Number of samples must be at least 1.")
  end
  return Vector[normalize( randn(dist.dim)+ 1im*randn(dist.dim) ) for _ in 1:n]
end

"""
Type corresponding to distribution of mixed states obtained by
tracing out part of a pure state obtained from a FubiniStudy
distribution.
"""
type FubiniStudyMixedState
    dim::Int64
    bath_dim::Int64
end

function rand(dist::FubiniStudyMixedState)
  return trace(projector(rand(FubiniStudyPureState(dist.dim+dist.bath_dim))),[dist.dim,dist.bath_dim],2)
end

"""
Type corresponding to the unitarily invariant distribution of
unitary transformations for some Hilbert space dimension.
"""
type HaarUnitary
    dim::Int64
end

function rand(dist::HaarUnitary)
  return svd(randn(dist.dim,dist.dim)+1im*randn(dist.dim,dist.dim))[1]
end

"""
Type corresponding to the unitarily invariant distribution of
unitary transformations for system and bath, such that the bath
(initially in the ground state) is traced out.
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
  ru = rand(HaarUnitary(dist.dim+dist.bath_dim)) * kron(eye(dist.dim),ket(0,dist.bath_dim))
  k = Matrix{Complex128}[]
  for ii=1:de
    push!(k, kron(eye(d),bra(ii-1,de)) * ru )
  end
  kraus2liou(k)
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
    rh = rand(Wigner(2.0),dist.dim)
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
