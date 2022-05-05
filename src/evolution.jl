"""
```
struct EvolutionSystem{Op}
    eval_tstep::Op
    grid::ExtendableGrid
    jac::SparseMatrixCSC{Float64, Int64}
    colors::Vector{Int}
end
```

Evolution system structure
"""
struct EvolutionSystem{Op}
    eval_tstep::Op
    grid::ExtendableGrid
    jac::SparseMatrixCSC{Float64, Int64}
    colors::Vector{Int}
end

"""
    EvolutionSystem(grid,implicit_tstep; jac=nothing)

Create evolution system on grid.
The signature of the callback is
    
    eval_tstep(residual,u,uold,grid,tstep)
"""
function EvolutionSystem(grid,eval_tstep; jac=nothing)
    N=num_nodes(grid)
    uold=rand(N)
    if isnothing(jac)
        A_h!(y,u)=eval_tstep(y,u,uold,grid,0.1)
        input = ones(N)
        output = similar(input)
        @time sparsity_pattern=Symbolics.jacobian_sparsity(A_h!,output,input)
        jac=Float64.(sparsity_pattern)
    end
    colors = SparseDiffTools.matrix_colors(jac)
    EvolutionSystem(eval_tstep,grid,jac,colors)
end




"""
     stdsparse(grid)
Create sparse matrix corresponding to P1 FEM for scalar problem.
"""
function stdsparse(grid)
    cn=grid[CellNodes]
    nn=num_nodes(grid)
    nc=num_cells(grid)
    m=ExtendableSparseMatrix(nn,nn)
    for ic=1:nc
        m[cn[1,ic],cn[1,ic]]=1
        m[cn[2,ic],cn[2,ic]]=1
        m[cn[3,ic],cn[3,ic]]=1

        m[cn[1,ic],cn[2,ic]]=1
        m[cn[2,ic],cn[1,ic]]=1
        m[cn[1,ic],cn[3,ic]]=1
        m[cn[3,ic],cn[1,ic]]=1
        m[cn[2,ic],cn[3,ic]]=1
        m[cn[3,ic],cn[2,ic]]=1
    end
    flush!(m)
    m.cscmatrix
end

"""      evolve(system::EvolutionSystem,inival;tend=1.0,nsteps=10)

Perform n timesteps in interval `(0,tend)`
"""
function evolve(system::EvolutionSystem,inival;tend=1.0,nsteps=10)
    t=0.0
    tsol=VoronoiFVM.TransientSolution(t,inival)
    Δt=tend/nsteps
    f!(y,u)=system.eval_tstep(y,u,inival,system.grid,Δt)
    j!(jac,u)=SparseDiffTools.forwarddiff_color_jacobian!(jac,f!,u,colorvec = system.colors)	
    df = OnceDifferentiable(f!, j!, copy(inival), copy(inival), system.jac)
    for istep=1:nsteps
	t=t+Δt
	fill!(system.jac,0)
        res=nlsolve(df,inival,method=:newton)
	append!(tsol,t,res.zero)
	inival=res.zero
    end
    tsol.t[end]=tend
    tsol
end

