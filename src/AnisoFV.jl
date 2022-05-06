module AnisoFV
using ExtendableGrids, ExtendableSparse, Triangulate
using LinearAlgebra,SparseArrays
using Symbolics,SparseDiffTools
using Tensors
using VoronoiFVM
using NLsolve
using Polynomials
using StaticArrays
    
include("fvmlib.jl")
export fvnorms, fvm_assemble!




include("femlib.jl")
export fenorms, fem_assemble!, fem_solve


include("afvmlib.jl")
export femfactors, baryfactors, circumfactors, afvm_solve

include("evolution.jl")
export EvolutionSystem,evolve,stdsparse,statsolve

include("operators.jl")
export ∇Λ∇, finitebell,ΛMatrix


end # module

