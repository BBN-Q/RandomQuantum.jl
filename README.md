# RandomQuantum

This package provides functions to sample from various random matrix ensebles
associated with quantum information applications.

Here is a table of currently support ensembles:

Type    | Common name
--------|----------
`GinibreEnsemble` | Ginibre Unitary Ensemble
`FubiniStudyPureState` | Fubini-Study ensemble
`FubiniStudyMixedState` | none
`HilbertSchmidtMixedState` | none
`BuresMixedState` | Bures Ensemble
`HaarUnitary` | Circular Unitary Ensemble (CUE) 
`HaarCPTPMap` | none
`GUEUnitary` | none
`GUECPTPMap` | none

The interface is emulates the interface of `Distributions.jl`,
although there is a lot missing at the moment.

# Examples

```julia
julia> using RandomQuantum

julia> rand(BuresMixedState(2))
2x2 Array{Complex{Float64},2}:
   0.771511+0.0im       -0.0632581+0.116198im
 -0.0632581-0.116198im    0.228489+0.0im    
```

# License

MIT License

# Copyright

(c) Raytheon BBN Technologies, 2015

# Contributors

Marcus P S (@marcusps)