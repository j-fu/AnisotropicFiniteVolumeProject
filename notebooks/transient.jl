### A Pluto.jl notebook ###
# v0.19.4

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
    using Revise
	Pkg.activate(joinpath(@__DIR__,".."))
	using AnisoFV
    using PlutoUI
    EL=PlutoUI.ExperimentalLayout
	using GridVisualize
	using ExtendableGrids
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

# ╔═╡ 411fa012-fed8-4aec-ad64-544753b75f95
md"""
## Parameters
"""

# ╔═╡ 80f88818-a342-4f5d-8a0e-1be36363a268
h=0.5

# ╔═╡ 329f3e2d-0105-482c-997a-79a69a7980f1
ugrid=usquare(;h)

# ╔═╡ 960f5bce-4a24-419d-8e2e-ec3c89861feb
rgrid=rsquare(;h)

# ╔═╡ 0b5709b9-3876-4967-b3f7-111f4ec1b31b
Λ=ΛMatrix(10,π/6)

# ╔═╡ 944e6bc5-6192-468e-b702-66e9da054a3b
md"""
## Scheme 1
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

# ╔═╡ 4855ce18-5c8d-419c-9e8a-942e26ab7250
tsol=evolve(sys,map(finitebell,grid);tend,nsteps=200);

# ╔═╡ 3034215b-e069-4421-942a-f572785b56fd
vis=GridVisualizer(dim=2,size=(300,300),colormap=:summer,levels=-0.1:0.05:1)

# ╔═╡ 3869679c-2011-404b-bec4-725afb76443b
scalarplot!(vis,grid,tsol(t),show=true)

# ╔═╡ c9c9a576-4713-42d4-b654-23f827696031
function minplot(tsol)
	mins=[minimum(tsol[i]) for i=1:length(tsol.t)]
	scalarplot(tsol.t,mins,size=(600,200))
end


# ╔═╡ 8d4c1451-160a-4932-ab0b-1cd49e086c7e
minplot(tsol)

# ╔═╡ b03775f9-d149-455a-a7a3-99acc46c8db7
md"""
## Scheme 2
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
vis2=GridVisualizer(dim=2,size=(300,300),colormap=:summer,levels=-0.1:0.1:1)

# ╔═╡ 1b91cad4-0abb-4ae5-bd78-dc05b38d6c04
scalarplot!(vis2,grid,tsol2(t2),show=true)

# ╔═╡ b977e12b-b709-4602-82b4-1b3103acbc75
minplot(tsol2)

# ╔═╡ 9f4321f5-302d-4236-82ae-95e0937e11ac
md"""
## Scheme 3
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
vis3=GridVisualizer(dim=2,size=(300,300),colormap=:summer,levels=-0.1:0.1:1)

# ╔═╡ 5ea03814-2c7f-4a4f-94e8-f32b2b322e63
scalarplot!(vis3,grid,tsol3(t3),show=true)

# ╔═╡ cfef3b64-1473-4cfb-84f1-ba2e3eb23eec
minplot(tsol3)

# ╔═╡ e1c83f34-6c42-49ba-b803-56ec9f069060
md"""
## Scheme 4
"""

# ╔═╡ e40f4140-9750-446b-9f9f-12f7b05b9754
xsqrt(x)= x<0 ? 0 : sqrt(x+1.0e-20)

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
vis4=GridVisualizer(dim=2,size=(300,300),colormap=:summer,levels=-0.1:0.1:1)

# ╔═╡ da0a0807-39f3-4dec-9756-18ced2aad16d
scalarplot!(vis4,grid,tsol4(t4),show=true)

# ╔═╡ 4a77b190-8d6b-4204-a8ab-1e3f7f66769e
minplot(tsol4)

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═5a4c1695-28be-428a-9983-f94b3fa38160
# ╠═7bda3f49-2d03-4e81-ac8a-b4f80db0ef5d
# ╟─411fa012-fed8-4aec-ad64-544753b75f95
# ╠═80f88818-a342-4f5d-8a0e-1be36363a268
# ╠═329f3e2d-0105-482c-997a-79a69a7980f1
# ╠═960f5bce-4a24-419d-8e2e-ec3c89861feb
# ╠═0b5709b9-3876-4967-b3f7-111f4ec1b31b
# ╟─944e6bc5-6192-468e-b702-66e9da054a3b
# ╠═2938b43a-aefd-402c-beb4-17cfe8b6c870
# ╠═639724fa-a638-4cfd-a632-497ac8a901cf
# ╠═aef21e27-4964-45f3-b215-972599114763
# ╠═00a4869f-57cb-4772-a741-d8016eb471d2
# ╠═9f1da51b-7076-4043-8906-276c1ad76fa5
# ╠═4855ce18-5c8d-419c-9e8a-942e26ab7250
# ╠═3034215b-e069-4421-942a-f572785b56fd
# ╠═3869679c-2011-404b-bec4-725afb76443b
# ╠═8d4c1451-160a-4932-ab0b-1cd49e086c7e
# ╠═c9c9a576-4713-42d4-b654-23f827696031
# ╟─b03775f9-d149-455a-a7a3-99acc46c8db7
# ╠═de9145dd-9747-40a0-b00e-9f2c94286a2e
# ╠═7795f2b5-1609-4950-a269-a50608b82496
# ╠═aece7bbb-4944-48cd-b8e4-aaf207f68067
# ╠═1d77b66b-748d-423f-a13d-b601fb8cc735
# ╠═280fb497-dbba-43d5-b2e4-1e865f82dae2
# ╠═6f508fc6-9801-4340-89f9-3ee12e6c1e05
# ╠═1b91cad4-0abb-4ae5-bd78-dc05b38d6c04
# ╠═b977e12b-b709-4602-82b4-1b3103acbc75
# ╟─9f4321f5-302d-4236-82ae-95e0937e11ac
# ╠═29bbf6c4-ad4e-4024-808e-2cf2344533af
# ╠═1c4a3d38-2adb-4b5a-b5f5-a944fb1da93a
# ╠═e3aabc6f-f279-4ff1-ad49-bd2e0066f188
# ╠═3a329b91-1592-4b5d-92a4-77e8532e154c
# ╠═790c47a7-650e-4ced-a056-a96a37f30993
# ╠═5ea03814-2c7f-4a4f-94e8-f32b2b322e63
# ╠═cfef3b64-1473-4cfb-84f1-ba2e3eb23eec
# ╟─e1c83f34-6c42-49ba-b803-56ec9f069060
# ╠═e40f4140-9750-446b-9f9f-12f7b05b9754
# ╠═a42a8dfe-9bec-4fcd-8ba0-6b64f67add88
# ╠═34281644-f132-4b63-91aa-1a67d2487213
# ╠═081d5177-df1e-4530-ab11-db9940126717
# ╠═c15d2afc-4e54-48f4-9288-d4f482ce3aa3
# ╠═54afbd40-6ad0-4071-9abf-c5c1443f6e25
# ╠═da0a0807-39f3-4dec-9756-18ced2aad16d
# ╠═4a77b190-8d6b-4204-a8ab-1e3f7f66769e
