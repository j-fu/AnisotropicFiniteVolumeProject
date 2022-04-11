using Pkg
Pkg.activate(joinpath(@__DIR__,".."))

using VoronoiFVM
using ExtendableGrids
using GridVisualize
using PyPlot


## Flux function which describes the flux
## between neigboring control volumes
function laplace_flux!(f,u,edge)
    f[1]=u[1,1]-u[1,2]
end


function laplace_main(;n=20)
    nspecies=1 
    ispec=1    
    X=collect(0:1.0/n:1)
    grid=VoronoiFVM.Grid(X,X)


    physics=VoronoiFVM.Physics(flux=laplace_flux!)
    sys=VoronoiFVM.System(grid,physics)
    enable_species!(sys,ispec,[1])
    boundary_dirichlet!(sys,ispec,1,0.0)
    boundary_dirichlet!(sys,ispec,3,1.0)
    inival=unknowns(sys,inival=0)
    solution=unknowns(sys)
    solve!(solution,inival,sys)
    nf=nodeflux(sys,solution)
    scalarplot(grid,solution[1,:],Plotter=PyPlot)
end

