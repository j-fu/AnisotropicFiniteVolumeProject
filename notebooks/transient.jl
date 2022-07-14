### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	using Pkg
        ENV["LC_NUMERIC"]="C"
	Pkg.activate(joinpath(@__DIR__,".."))
	using AnisotropicFiniteVolumeProject
    using Revise
    using PlutoUI
	using PyPlot
    EL=PlutoUI.ExperimentalLayout
	using GridVisualize
	using ExtendableGrids
	using DrWatson
	using SimplexGridFactory,Triangulate
	import PlutoVista
	default_plotter!(PlutoVista)
	using LinearAlgebra
	TableOfContents()
end

# ╔═╡ 5a4c1695-28be-428a-9983-f94b3fa38160
function usquare(;h=0.1)
    # Create a SimplexGridBuilder structure which can collect
    # geometry information
    builder=SimplexGridBuilder(Generator=Triangulate)
    scale=10
    # Add points, record their numbers
    p1=point!(builder,-1*scale,-1*scale)
    p2=point!(builder,1*scale,-1*scale)
    p3=point!(builder,1*scale,1*scale)
    p4=point!(builder,-1*scale,1*scale)
    
    # Connect points by respective facets (segments)
    facetregion!(builder,1)
    facet!(builder,p1,p2)
    facetregion!(builder,2)
    facet!(builder,p2,p3)
    facetregion!(builder,3)
    facet!(builder,p3,p4)
    facetregion!(builder,4)
    facet!(builder,p4,p1)
    options!(builder,maxvolume=0.1)
    simplexgrid(builder,maxvolume=0.75*h^2)
end

# ╔═╡ 7bda3f49-2d03-4e81-ac8a-b4f80db0ef5d
function rsquare(;h=0.1)
	X=-10:h:10
	simplexgrid(X,X)
end

# ╔═╡ ea1e02a5-e4a7-445e-aa5c-4cc2dd0b51a4
plottimes=[0.1,0.3,0.6,1.0]

# ╔═╡ 32cec990-461a-45f1-8e34-b3f0b9451d6d
function mysave(fname,v)
    GridVisualize.save(plotsdir(fname),v)
	@info """saved to $(plotsdir(fname))"""
end

# ╔═╡ c825f9d6-e040-4a5e-87de-4de3282e7023
function fourplots(grid,tsol; fname="test.png", times=plottimes)
	vis=GridVisualizer(layout=(2,2),resolution=(700,600),Plotter=PyPlot)
	myplot(i,j,t)=
	scalarplot!(vis[i,j],grid,tsol(t),colormap=:summer,levels=5, title="t=$t")	
	myplot(1,1,times[1])
	myplot(1,2,times[2])
	myplot(2,1,times[3])
	myplot(2,2,times[4])
		
	mysave(fname,vis)
	reveal(vis)
end

# ╔═╡ 411fa012-fed8-4aec-ad64-544753b75f95
md"""
## Parameters
"""

# ╔═╡ 80f88818-a342-4f5d-8a0e-1be36363a268
h=0.3

# ╔═╡ ee3fa050-a2f5-4b53-97b2-0fd9c95b8935
md"""
Create a completely unstructured grid:
"""

# ╔═╡ 329f3e2d-0105-482c-997a-79a69a7980f1
ugrid=usquare(;h)

# ╔═╡ 2c380b32-d3b4-415c-9e81-4d105391ea27
md"""
Create a triangular grid from structured grid;
"""

# ╔═╡ 960f5bce-4a24-419d-8e2e-ec3c89861feb
rgrid=rsquare(;h)

# ╔═╡ 0b5709b9-3876-4967-b3f7-111f4ec1b31b
Λ=ΛMatrix(100,π/6)

# ╔═╡ c9c9a576-4713-42d4-b654-23f827696031
function minplot(tsol)
	mins=[minimum(tsol[i]) for i=1:length(tsol.t)]
	scalarplot(tsol.t,mins,size=(600,200))
end

# ╔═╡ 944e6bc5-6192-468e-b702-66e9da054a3b
md"""
## Scheme 1
"""

# ╔═╡ d03f3d5d-70d6-4cd0-8fa5-735844b2d454
md"""
This schem corresponds to the simple linear finite volume scheme for the equation

```math
\partial_t u  - \nabla\cdot \Lambda  \nabla u = 0
```
"""

# ╔═╡ 2938b43a-aefd-402c-beb4-17cfe8b6c870
function diffusion_step(f,u,uold,sys,Δt)
	e=sys.e
	ω=sys.ω
	tris::Matrix{Int64}=sys.grid[CellNodes]
	ntri=size(e,2)
	en=[2 3 ; 3 1 ; 1 2 ]'
	fill!(f,0.0)
	for itri=1:ntri
		for iedge=1:3
			k=tris[en[1,iedge],itri]		
			l=tris[en[2,iedge],itri]	
			du=e[iedge,itri]*(u[k]-u[l])
			f[k]+=du
			f[l]-=du
		end
	for inode=1:3
			k=tris[inode,itri]
			f[k]+=ω[inode,itri]*(u[k]-uold[k])/Δt
		end
	end
end

# ╔═╡ 639724fa-a638-4cfd-a632-497ac8a901cf
grid=ugrid

# ╔═╡ aef21e27-4964-45f3-b215-972599114763
sys=EvolutionSystem(grid,diffusion_step;jac=stdsparse(grid),Λ);

# ╔═╡ 00a4869f-57cb-4772-a741-d8016eb471d2
tend=10

# ╔═╡ 9f1da51b-7076-4043-8906-276c1ad76fa5
@bind t Slider(0:0.01:tend,show_value=true)

# ╔═╡ b7ae126b-651d-4a25-8096-be96b16a6f0e
md"""
We calculate the evolution of our finitebell function
"""

# ╔═╡ 4855ce18-5c8d-419c-9e8a-942e26ab7250
tsol=evolve(sys,map(finitebell,grid);tend,nsteps=200);

# ╔═╡ 3034215b-e069-4421-942a-f572785b56fd
vis=GridVisualizer(dim=2,size=(300,300),colormap=:summer,levels=5)

# ╔═╡ 3869679c-2011-404b-bec4-725afb76443b
scalarplot!(vis,grid,tsol(t),show=true)

# ╔═╡ 8d4c1451-160a-4932-ab0b-1cd49e086c7e
minplot(tsol)

# ╔═╡ 0e06aa3a-9a5b-4744-b3eb-22df1492c9a8
minimum(tsol)

# ╔═╡ b2fc28d2-7dd7-4246-b969-4bf9fe7e0df9
fourplots(grid,tsol,fname="transient1.png")

# ╔═╡ c340af0e-8a29-4fc8-81c9-d0451367250a
md"""
As we see for scheme 1,  though we should expect a nonnegative solution due the the nonnegative initial value, we find negative solution values at the begining of the evolution.
"""

# ╔═╡ b03775f9-d149-455a-a7a3-99acc46c8db7
md"""
## Scheme 2
"""

# ╔═╡ b219629b-9c97-4598-8848-3809854d0752
md"""
Here, we solve 
```math
\partial_t u  - \nabla \cdot \Lambda u \nabla \log u = 0
```
using maximum/minimum of  `u` depending on the sign of the entries in the local stiffness matrix.
"""

# ╔═╡ de9145dd-9747-40a0-b00e-9f2c94286a2e
xlog(u)= u<1.0e-20 ? log(1.0e-20) : log(u)

# ╔═╡ 7795f2b5-1609-4950-a269-a50608b82496
function diffusion_step2(f,u,uold,sys,Δt)
	e=sys.e
	ω=sys.ω
	tris::Matrix{Int64}=sys.grid[CellNodes]
	ntri=size(e,2)
	en=[2 3 ; 3 1 ; 1 2 ]'
	fill!(f,0.0)
	for itri=1:ntri
		for iedge=1:3
			k=tris[en[1,iedge],itri]		
			l=tris[en[2,iedge],itri]	
			ukl=e[iedge,itri]>0  ? max(u[k],u[l]) : min(u[k],u[l])
			du=ukl*e[iedge,itri]*(xlog(u[k])-xlog(u[l]))
			f[k]+=du
			f[l]-=du
		end
	for inode=1:3
			k=tris[inode,itri]
			f[k]+=ω[inode,itri]*(u[k]-uold[k])/Δt
		end
	end
end

# ╔═╡ aece7bbb-4944-48cd-b8e4-aaf207f68067
sys2=EvolutionSystem(grid,diffusion_step2;jac=stdsparse(grid),Λ);

# ╔═╡ 1d77b66b-748d-423f-a13d-b601fb8cc735
tsol2=evolve(sys2,map(finitebell,grid);tend,nsteps=200);

# ╔═╡ 280fb497-dbba-43d5-b2e4-1e865f82dae2
@bind t2 Slider(0:0.01:tend,show_value=true)

# ╔═╡ 6f508fc6-9801-4340-89f9-3ee12e6c1e05
vis2=GridVisualizer(dim=2,size=(300,300),colormap=:summer,levels=5)

# ╔═╡ 1b91cad4-0abb-4ae5-bd78-dc05b38d6c04
scalarplot!(vis2,grid,tsol2(t2),show=true)

# ╔═╡ b977e12b-b709-4602-82b4-1b3103acbc75
minplot(tsol2)

# ╔═╡ cdc0c1e1-a72d-4047-8a9d-04ad1c984d39
minimum(tsol2)

# ╔═╡ 1b57391b-6abb-4c8e-b7fe-c83b7859c272
md"""
As we see, the solution stays positive.
"""

# ╔═╡ b48b4d86-7a17-45a5-8153-de58e86c21c8
fourplots(grid,tsol2,fname="transient2.png")

# ╔═╡ 9f4321f5-302d-4236-82ae-95e0937e11ac
md"""
## Scheme 3
"""

# ╔═╡ d58ff61e-44a9-4906-8808-29dae7090ea0
md"""
Here we solve the same problem as with scheme2, but use an upwinding
depnding on the flux 
"""

# ╔═╡ 29bbf6c4-ad4e-4024-808e-2cf2344533af
function diffusion_step3(f,u,uold,sys,Δt)
	e=sys.e
	ω=sys.ω
	tris::Matrix{Int64}=sys.grid[CellNodes]
	ntri=size(e,2)
	en=[2 3 ; 3 1 ; 1 2 ]'
	fill!(f,0.0)
	for itri=1:ntri
		for iedge=1:3
			k=tris[en[1,iedge],itri]		
			l=tris[en[2,iedge],itri]
			Fkl=e[iedge,itri]*(xlog(u[k])-xlog(u[l]))
			ukl=Fkl>0  ? u[k] : u[l]
			du=ukl*Fkl
			f[k]+=du
			f[l]-=du
		end
	for inode=1:3
			k=tris[inode,itri]
			f[k]+=ω[inode,itri]*(u[k]-uold[k])/Δt
		end
	end
end

# ╔═╡ 1c4a3d38-2adb-4b5a-b5f5-a944fb1da93a
sys3=EvolutionSystem(grid,diffusion_step3;jac=stdsparse(grid),Λ);

# ╔═╡ e3aabc6f-f279-4ff1-ad49-bd2e0066f188
tsol3=evolve(sys3,map(finitebell,grid);tend,nsteps=200);

# ╔═╡ 3a329b91-1592-4b5d-92a4-77e8532e154c
@bind t3 Slider(0:0.01:tend,show_value=true)

# ╔═╡ 790c47a7-650e-4ced-a056-a96a37f30993
vis3=GridVisualizer(dim=2,size=(300,300),colormap=:summer,levels=5)

# ╔═╡ 5ea03814-2c7f-4a4f-94e8-f32b2b322e63
scalarplot!(vis3,grid,tsol3(t3),show=true)

# ╔═╡ cfef3b64-1473-4cfb-84f1-ba2e3eb23eec
minplot(tsol3)

# ╔═╡ 8a510441-81cc-4c41-8f21-2adda9fec5f1
minimum(tsol3)

# ╔═╡ 256bf385-5cdf-437f-b41b-91b4ddc54925
fourplots(grid,tsol3,fname="transient3.png")

# ╔═╡ 450a3c5f-93d9-43d6-9358-e05c5de4f76e
md"""
Once again, we see the positivity of the flux
"""

# ╔═╡ e1c83f34-6c42-49ba-b803-56ec9f069060
md"""
## Scheme 4
(incorrect version of Quenjel scheme, we skip this)
"""

# ╔═╡ e40f4140-9750-446b-9f9f-12f7b05b9754
xsqrt(x)= x<0 ? 0 : sqrt(x+1.0e-20)

# ╔═╡ 5c2d4fd2-1de8-4727-9a78-c3325d6565c7
md"""
Here, we use the same splitting ``\Lambda u = \sqrt(u)\Lambda\sqrt(u)`` as in the last scheme in the stationary case. In the stationary case, this has the best convergence rate of the upwinded schemes handling anisitropy. Here, we check it for nonnegativity.
"""

# ╔═╡ f29e3368-b59f-4e98-ab17-99f98af3bd13


# ╔═╡ a42a8dfe-9bec-4fcd-8ba0-6b64f67add88
function diffusion_step4(f,u,uold,sys,Δt)
	tris::Matrix{Int64}=sys.grid[CellNodes]
	coord::Matrix{Float64}=sys.grid[Coordinates]
	ntri=size(tris,2)
	en=[2 3 ; 3 1 ; 1 2 ]'
	fill!(f,0.0)
	ϕij=ones(eltype(u),3)
	ϕi=ones(eltype(u),3)
    for itri=1:ntri
		for inode=1:3
			k=tris[inode,itri]
			ϕi[inode]=xsqrt(u[k])
		end
		for iedge=1:3
			i=en[1,iedge]		
			j=en[2,iedge]
			ϕij[iedge]=0.5*(ϕi[i]+ϕi[j])		
		end
		ω,e=baryfactors(itri,ϕij,Λ,coord,tris)
		for iedge=1:3
			i=en[1,iedge]		
			j=en[2,iedge]
			k=tris[i,itri]		
			l=tris[j,itri]
			Fkl=e[iedge]*(xlog(u[k])-xlog(u[l]))
			ukl=Fkl>0  ? ϕi[i] : ϕi[j]
			du=ukl*Fkl
			f[k]+=du
			f[l]-=du
		end
	for inode=1:3
			k=tris[inode,itri]
			f[k]+=ω[inode]*(u[k]-uold[k])/Δt
		end
	end
end

# ╔═╡ 34281644-f132-4b63-91aa-1a67d2487213
sys4=EvolutionSystem(grid,diffusion_step4;jac=stdsparse(grid),Λ);

# ╔═╡ 081d5177-df1e-4530-ab11-db9940126717
tsol4=evolve(sys4,map(finitebell,grid);tend,nsteps=20);

# ╔═╡ c15d2afc-4e54-48f4-9288-d4f482ce3aa3
@bind t4 Slider(0:0.01:tend,show_value=true)

# ╔═╡ 54afbd40-6ad0-4071-9abf-c5c1443f6e25
vis4=GridVisualizer(dim=2,size=(300,300))

# ╔═╡ da0a0807-39f3-4dec-9756-18ced2aad16d
scalarplot!(vis4,grid,tsol4(t4),show=true,colormap=:summer,levels=5)

# ╔═╡ 4a77b190-8d6b-4204-a8ab-1e3f7f66769e
minplot(tsol4)

# ╔═╡ 4e5095a4-bef1-4d7b-b843-63b3bc956e99
minimum(tsol4)

# ╔═╡ da8fdc58-e534-4a55-9e34-694ece4ecb21
fourplots(grid,tsol4,fname="transient4.png")

# ╔═╡ 351923d7-da16-49d9-9b8a-e58e53dce044
md"""
As we see, this scheme delivers a nonnegative solution.
"""

# ╔═╡ 5b2af241-9d66-4b4f-841a-141c420958b1
md"""
## Scheme 5

This corresponds to the scheme in the paper of Quenjel
"""

# ╔═╡ 115c9d1e-fbd7-4c00-a23b-53106ea115f8
function diffusion_step5(f,u,uold,sys,Δt)
	tris::Matrix{Int64}=sys.grid[CellNodes]
	coord::Matrix{Float64}=sys.grid[Coordinates]
	ntri=size(tris,2)
	en=[2 3 ; 3 1 ; 1 2 ]'
	fill!(f,0.0)
	ϕij=ones(eltype(u),3)
	ϕi=ones(eltype(u),3)
    for itri=1:ntri
		for inode=1:3
			k=tris[inode,itri]
			ϕi[inode]=xsqrt(u[k])
		end
		for iedge=1:3
			i=en[1,iedge]		
			j=en[2,iedge]
			ϕij[iedge]=0.5*(ϕi[i]+ϕi[j])		
		end
        ω,a=baryfactorsx(itri,Λ,coord,tris)
		for iedge=1:3
			kl=iedge
			i=en[1,iedge]		
			j=en[2,iedge]
			k=tris[i,itri]		
			l=tris[j,itri]
			Fkl=0.0
			for jedge=1:3
				ii=tris[en[1,jedge],itri]
				iτ=tris[en[2,jedge],itri]
			    Fkl+=a[kl,jedge]*ϕij[jedge]*(xlog(u[ii])-xlog(u[iτ]))
			end			
			ukl=Fkl>0  ? ϕi[i] : ϕi[j]
			du=ukl*Fkl
			f[k]+=du
			f[l]-=du
		end
		
	for inode=1:3
			k=tris[inode,itri]
			f[k]+=ω[inode]*(u[k]-uold[k])/Δt
		end
	end
end

# ╔═╡ 89f6a669-16ec-4ac8-9f89-bfa322523fe4
sys5=EvolutionSystem(grid,diffusion_step5;jac=stdsparse(grid),Λ);

# ╔═╡ 7229f21e-8b92-4e8e-a67f-4c8c9c2c0aa7
tsol5=evolve(sys5,map(finitebell,grid);tend,nsteps=20);

# ╔═╡ 4a9d88c9-4492-401e-bf5f-43e24874cc19
@bind t5 Slider(0:0.01:tend,show_value=true)

# ╔═╡ 95ed710f-b0b2-43fc-ab96-232166846972
vis5=GridVisualizer(dim=2,size=(300,300))

# ╔═╡ 7702b595-bf66-4913-a584-cab5c77ebaa4
scalarplot!(vis5,grid,tsol5(t5),show=true,colormap=:summer,levels=5)

# ╔═╡ 454ee521-6581-425d-81b3-8d0a116f2e87
minplot(tsol5)

# ╔═╡ 01234ee8-c6c0-45c9-9bab-d306940308d8
minimum(tsol5)

# ╔═╡ 68106561-3e5a-40db-8516-f3acd4383122
fourplots(grid,tsol5,fname="transient5.png")

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═5a4c1695-28be-428a-9983-f94b3fa38160
# ╠═7bda3f49-2d03-4e81-ac8a-b4f80db0ef5d
# ╠═ea1e02a5-e4a7-445e-aa5c-4cc2dd0b51a4
# ╠═32cec990-461a-45f1-8e34-b3f0b9451d6d
# ╠═c825f9d6-e040-4a5e-87de-4de3282e7023
# ╟─411fa012-fed8-4aec-ad64-544753b75f95
# ╠═80f88818-a342-4f5d-8a0e-1be36363a268
# ╟─ee3fa050-a2f5-4b53-97b2-0fd9c95b8935
# ╠═329f3e2d-0105-482c-997a-79a69a7980f1
# ╟─2c380b32-d3b4-415c-9e81-4d105391ea27
# ╠═960f5bce-4a24-419d-8e2e-ec3c89861feb
# ╠═0b5709b9-3876-4967-b3f7-111f4ec1b31b
# ╠═c9c9a576-4713-42d4-b654-23f827696031
# ╟─944e6bc5-6192-468e-b702-66e9da054a3b
# ╟─d03f3d5d-70d6-4cd0-8fa5-735844b2d454
# ╠═2938b43a-aefd-402c-beb4-17cfe8b6c870
# ╠═639724fa-a638-4cfd-a632-497ac8a901cf
# ╠═aef21e27-4964-45f3-b215-972599114763
# ╠═00a4869f-57cb-4772-a741-d8016eb471d2
# ╠═9f1da51b-7076-4043-8906-276c1ad76fa5
# ╟─b7ae126b-651d-4a25-8096-be96b16a6f0e
# ╠═4855ce18-5c8d-419c-9e8a-942e26ab7250
# ╠═3034215b-e069-4421-942a-f572785b56fd
# ╠═3869679c-2011-404b-bec4-725afb76443b
# ╠═8d4c1451-160a-4932-ab0b-1cd49e086c7e
# ╠═0e06aa3a-9a5b-4744-b3eb-22df1492c9a8
# ╠═b2fc28d2-7dd7-4246-b969-4bf9fe7e0df9
# ╟─c340af0e-8a29-4fc8-81c9-d0451367250a
# ╟─b03775f9-d149-455a-a7a3-99acc46c8db7
# ╟─b219629b-9c97-4598-8848-3809854d0752
# ╠═de9145dd-9747-40a0-b00e-9f2c94286a2e
# ╠═7795f2b5-1609-4950-a269-a50608b82496
# ╠═aece7bbb-4944-48cd-b8e4-aaf207f68067
# ╠═1d77b66b-748d-423f-a13d-b601fb8cc735
# ╠═280fb497-dbba-43d5-b2e4-1e865f82dae2
# ╠═6f508fc6-9801-4340-89f9-3ee12e6c1e05
# ╠═1b91cad4-0abb-4ae5-bd78-dc05b38d6c04
# ╠═b977e12b-b709-4602-82b4-1b3103acbc75
# ╠═cdc0c1e1-a72d-4047-8a9d-04ad1c984d39
# ╟─1b57391b-6abb-4c8e-b7fe-c83b7859c272
# ╠═b48b4d86-7a17-45a5-8153-de58e86c21c8
# ╟─9f4321f5-302d-4236-82ae-95e0937e11ac
# ╟─d58ff61e-44a9-4906-8808-29dae7090ea0
# ╠═29bbf6c4-ad4e-4024-808e-2cf2344533af
# ╠═1c4a3d38-2adb-4b5a-b5f5-a944fb1da93a
# ╠═e3aabc6f-f279-4ff1-ad49-bd2e0066f188
# ╠═3a329b91-1592-4b5d-92a4-77e8532e154c
# ╠═790c47a7-650e-4ced-a056-a96a37f30993
# ╠═5ea03814-2c7f-4a4f-94e8-f32b2b322e63
# ╠═cfef3b64-1473-4cfb-84f1-ba2e3eb23eec
# ╠═8a510441-81cc-4c41-8f21-2adda9fec5f1
# ╠═256bf385-5cdf-437f-b41b-91b4ddc54925
# ╟─450a3c5f-93d9-43d6-9358-e05c5de4f76e
# ╟─e1c83f34-6c42-49ba-b803-56ec9f069060
# ╠═e40f4140-9750-446b-9f9f-12f7b05b9754
# ╟─5c2d4fd2-1de8-4727-9a78-c3325d6565c7
# ╠═f29e3368-b59f-4e98-ab17-99f98af3bd13
# ╠═a42a8dfe-9bec-4fcd-8ba0-6b64f67add88
# ╠═34281644-f132-4b63-91aa-1a67d2487213
# ╠═081d5177-df1e-4530-ab11-db9940126717
# ╠═c15d2afc-4e54-48f4-9288-d4f482ce3aa3
# ╠═54afbd40-6ad0-4071-9abf-c5c1443f6e25
# ╠═da0a0807-39f3-4dec-9756-18ced2aad16d
# ╠═4a77b190-8d6b-4204-a8ab-1e3f7f66769e
# ╠═4e5095a4-bef1-4d7b-b843-63b3bc956e99
# ╠═da8fdc58-e534-4a55-9e34-694ece4ecb21
# ╟─351923d7-da16-49d9-9b8a-e58e53dce044
# ╟─5b2af241-9d66-4b4f-841a-141c420958b1
# ╠═115c9d1e-fbd7-4c00-a23b-53106ea115f8
# ╠═89f6a669-16ec-4ac8-9f89-bfa322523fe4
# ╠═7229f21e-8b92-4e8e-a67f-4c8c9c2c0aa7
# ╠═4a9d88c9-4492-401e-bf5f-43e24874cc19
# ╠═95ed710f-b0b2-43fc-ab96-232166846972
# ╠═7702b595-bf66-4913-a584-cab5c77ebaa4
# ╠═454ee521-6581-425d-81b3-8d0a116f2e87
# ╠═01234ee8-c6c0-45c9-9bab-d306940308d8
# ╠═68106561-3e5a-40db-8516-f3acd4383122
