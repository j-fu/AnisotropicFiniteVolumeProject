module AnisoFV
using ExtendableGrids, ExtendableSparse
using LinearAlgebra,SparseArrays
using Symbolics,SparseDiffTools
using VoronoiFVM
using NLsolve


include("fvmlib.jl")
export fvnorms, fvm_assemble!



include("femlib.jl")
export fenorms, fem_assemble!


include("evolution.jl")
export EvolutionSystem,evolve,stdsparse




end # module

