### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
    using Revise
	using AnisotropicFiniteVolumeProject
    using Polynomials
	
    using PlutoUI
    EL=PlutoUI.ExperimentalLayout
	using GridVisualize
	using PyPlot
	using VoronoiFVM
	using ExtendableGrids
	using DrWatson
	default_plotter!(PyPlot)
	PyPlot.svg(true)
	using LinearAlgebra
	TableOfContents()
end

# ╔═╡ de8b10df-42d6-43a5-b9d7-59bf7a4d99f3
function mysave(fname,v)
    GridVisualize.save(plotsdir(fname),v)
	@info """saved to $(plotsdir(fname))"""
end

# ╔═╡ 223d2526-9e76-43c4-8ba6-78b56342b56c
function mysave(fname)
    savefig(plotsdir(fname))
	@info """saved to $(plotsdir(fname))"""
end

# ╔═╡ 9eb0b277-5a96-436e-a233-051edfa7cb24
md"""
# Stationary anisotropic Problem
"""

# ╔═╡ acaed9cb-5ed0-4137-b65c-904b26b4848f
md"""
## Tools to define a stationary test problem

### Twice differentiable function with finite support
"""

# ╔═╡ 207a355d-ff8f-4404-90db-daf9341f4f98
@doc finitebell_core

# ╔═╡ 33755343-8709-4c88-95f9-9758f1d61f26
@doc finitebell

# ╔═╡ e5a25fb2-4525-4ca4-b155-b00f90fea311
@doc d1finitebell

# ╔═╡ 6b417d68-7428-4ebe-8932-eec7424f100c
@doc d2finitebell

# ╔═╡ cfc41a2a-3df9-4fc2-82ec-2769b25e9d56
let 
	X=-1.5:0.01:1.5
	v=GridVisualizer(size=(600,300),legend=:lb)
    scalarplot!(v,X,finitebell.(X),color=:red,label="finitebell")
    scalarplot!(v,X,d1finitebell.(X),clear=false,color=:green,label="d1finitebell")
    scalarplot!(v,X,d2finitebell.(X),clear=false,color=:blue,label="d2finitebell")
	mysave("finitebell1d.png",v)
	reveal(v)
end

# ╔═╡ 093cada5-9eba-40e1-bff7-d7aee5cf6b19
md"""
### ∇Λ∇ operator
"""

# ╔═╡ d64a304c-e970-4c1a-b2d6-8314bd69df98
@doc ∇Λ∇

# ╔═╡ 9669a6e5-0593-49ed-9fa6-35a1e3986548
md"""
This operator can be used to calcualate the Laplacian of the finitebell function
which is visualized below.
"""

# ╔═╡ 344b51fb-2a23-4620-ac30-4129e7e48df7
let
	X=-1.5:0.01:1.5
	g2=simplexgrid(X,X)
	fb=map(finitebell, g2)
	d2fb=map((x,y)->∇Λ∇(finitebell,[x,y]), g2)

	vis=GridVisualizer(layout=(1,2),size=(600,300))
	scalarplot!(vis[1,1],g2,fb,title="finitebell",limits=(-1,1),colormap=:bwr)
	scalarplot!(vis[1,2],g2,d2fb,title="Δ finitebell",limits=(-10,10),colormap=:bwr)
		mysave("finitebell2d.png",vis)
	reveal(vis)
end

# ╔═╡ 04838dde-4924-4b67-b9b2-1ecd534fe861
@doc ΛMatrix

# ╔═╡ 86c83f88-bf48-48c9-ba57-b6dae3574c18
ΛMatrix(10,0)

# ╔═╡ 101e6962-1046-4f98-a363-a666c026b95e
ΛMatrix(10,π/2)

# ╔═╡ 06abe474-3862-440e-9e02-59fb897b1562
ΛMatrix(10,π/4)

# ╔═╡ 619675e5-abea-48e1-83ef-3f0ff4395eed
md"""
### Functions to create grids
"""

# ╔═╡ 54af9427-239b-4cc2-adaf-c443edb90df3
@doc rgrid

# ╔═╡ 250cbba5-334a-494c-9d03-f72e8fb38937
@doc tgrid

# ╔═╡ 3e191aec-1121-4584-99e9-953c894fb13c
let
	v=GridVisualizer(size=(600,300),layout=(1,2))
    gridplot!(v[1,1],rgrid(h=0.1,scale=1.1),linewidth=0.5,title="rgrid")
     gridplot!(v[1,2],tgrid(h=0.1,scale=1.1),linewidth=0.5,title="tgrid")
	mysave("grids.png",v)
	reveal(v)
end

# ╔═╡ 7b273b52-be4f-4a32-9561-d4abfb79761b
md"""
## Test problem

```math
    -\nabla \cdot (\Lambda \nabla u)=f
```
with homogeneous Dirichlet boundary condition. We choose `f` sucht that the exact solution is our `finitebell` function.
"""

# ╔═╡ 7738bc89-efd4-4b37-b5b4-0de400657f7e
exact(grid)=map((x,y)->finitebell(x,y),grid)

# ╔═╡ a20f6356-f05f-4927-a0dd-6226b33c1433
maxlev=5

# ╔═╡ ed8a0ecd-5563-4351-8291-b35cabfc1b03
function fourtests(test::Function,fname;h=0.1,scale=1.1)
	vis=GridVisualizer(layout=(2,2),resolution=(700,600))
	rg=rgrid(;h,scale)
	tg=tgrid(;h,scale)
	scalarplot!(vis[1,1],rg,test(rg,ΛMatrix(100,0)),colormap=:bwr,flimits=(-1,1),levels=-1:0.2:1,title="rgrid,ΛMatrix(100,0)")
	scalarplot!(vis[2,1],rg,test(rg,ΛMatrix(100,π/4)),colormap=:bwr,flimits=(-1,1),levels=-1:0.2:1,title="rgrid,ΛMatrix(100,π/4)")
	scalarplot!(vis[1,2],tg,test(tg,ΛMatrix(100,0)),colormap=:bwr,flimits=(-1,1),levels=-1:0.2:1,title="tgrid,ΛMatrix(100,0)")
	scalarplot!(vis[2,2],tg,test(tg,ΛMatrix(100,π/4)),colormap=:bwr,flimits=(-1,1),levels=-1:0.2:1,title="tgrid,ΛMatrix(100,π/4)")
	mysave(fname,vis)
	reveal(vis)
end

# ╔═╡ 687f9d1e-4dc1-4791-a4c4-28f7e6b46983
H=[0.2*2.0^-i for i=0:maxlev-1]

# ╔═╡ 50809264-3603-41d9-923c-5a163f221344
function fourconv(test::Function;scale=1.1,maxlev=3)
	
	xnorm(g,u)=fenorms(u,g[Coordinates],g[CellNodes])
	function gtest(xgrid,Λ,H)
		norms=[]
		for h∈H
			g=xgrid(;scale,h)
	  	    push!(norms,xnorm(g,test(g,Λ)-exact(g)))
		end
		norms
	end
    gtest(rgrid,ΛMatrix(100,0),H),		
    gtest(rgrid,ΛMatrix(100,π/4),H),		
    gtest(tgrid,ΛMatrix(100,0),H),		
    gtest(tgrid,ΛMatrix(100,π/4),H)		

end

# ╔═╡ 34d148a4-8994-4f7f-9774-40a63644bcee
md"""
### Voronoi FVM
We define a VoronoiFVM.jl system and solve it. Here we project the diffusion tensor onto the grid edge.
"""

# ╔═╡ d86eef8a-99e1-4fc4-ba9f-213a41b86598
function VoronoiFVMTest(grid,Λ)

	function flux(y,u,edge)
		normal=edge[:,1]-edge[:,2]
		normal=normal/norm(normal)
		λ=norm(Λ*normal)
		y[1]=λ*(u[1,1]-u[1,2])
    end

	function source(y,node)
		y[1]=-∇Λ∇(finitebell,[node[1],node[2]],Λ)	
	end
	
	function bcondition(y,u,bnode)
		boundary_dirichlet!(y,u,bnode,value=0)		
	end
    sys=VoronoiFVM.System(grid;flux,bcondition,source,species=[1])
	sol=solve(sys)
	sol[1,:]
end

# ╔═╡ 4b439328-18d9-42fc-9575-85594a342b62
fourtests(VoronoiFVMTest,"voronoifvm.png")

# ╔═╡ 77ad9955-09bd-4042-9aad-0f94fa904ff7
VoronoiFVMNorms=fourconv(VoronoiFVMTest;maxlev)

# ╔═╡ 8b86f35e-ff6b-482c-b150-796da7db68ee
md"""
### FEM
"""

# ╔═╡ 49b5bc2a-96df-483d-862a-e0ee22f5e0b0
@doc fem_solve

# ╔═╡ 84ebe218-4dd0-440a-91ce-f37c5d06a341
function FEMTest(grid,Λ)
	f(x,y)=-∇Λ∇(finitebell,[x,y],Λ)
	β(x,y)=0.0
	fem_solve(grid,Λ,f,β)
end

# ╔═╡ 251eeb13-2730-413f-896f-5738fda9e801
fourtests(FEMTest,"fem.png")

# ╔═╡ 80ca622c-a1ed-4216-b25c-428e591ca32a
FEMNorms=fourconv(FEMTest;maxlev)

# ╔═╡ 97dc31e8-a4e5-4294-b214-da3e53089cc8
md"""
### "Simple" Barycentric finite volumes
"""

# ╔═╡ 648a48f0-2af4-4b9e-80f5-1926653c02ee
@doc EvolutionSystem

# ╔═╡ 5f53d081-97ca-4763-9b43-5822b97dd6a2
finitebellx(x,y)=finitebell(x,y)+0.01*rand()

# ╔═╡ 7511e717-3ad2-4d3b-8445-b11b5a5bb3fc
md"""
### Nonlinear BaryFVM, simple upwind
"""

# ╔═╡ 038dd7cb-07eb-481c-9459-29afab26c9bb
md"""
Logarithm cut of at negative values
"""

# ╔═╡ 1cfc228b-fa93-49e7-98c0-0f8d8dc45f8d
xlog(u)= u<1.0e-20 ? log(1.0e-20) : log(u)

# ╔═╡ 6fd997bd-1641-4d4b-a327-bea39cea4ba0
function homogeneous_dirichlet(f,u,grid)
    bfacenodes=grid[BFaceNodes]
    nbfaces=size(bfacenodes,2)
    bfaceregions=grid[BFaceRegions]
    for ibface in 1:nbfaces
        for idim=1:2
            i1=bfacenodes[idim,ibface];
            f[i1]+=1.0e30*u[i1]
        end
    end
end

# ╔═╡ 823f8838-a6ee-4ccc-8f10-27510499e01b
function BaryFVMTest1(grid,Λ)
# Λ is part of edge factors e!
function diffusion1(f,u,uold,sys,Δt)
	e=sys.e
	ω=sys.ω
	tris::Matrix{Int64}=sys.grid[CellNodes]
	coord::Matrix{Float64}=sys.grid[Coordinates]
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
			f[k]+=ω[inode,itri]*∇Λ∇(finitebell,coord[:,k],Λ)
		end
	end
 homogeneous_dirichlet(f,u,sys.grid)	
end
	
    sys=EvolutionSystem(grid,diffusion1;jac=stdsparse(grid),Λ=Λ);
	
	statsolve(sys,map(finitebellx,grid))
	
end

# ╔═╡ d5d60909-cf40-4b0f-9626-98e2b1fb853a
fourtests(BaryFVMTest1,"baryfvm1.png")

# ╔═╡ e758fc09-aa0a-4ed3-852a-faa30bc9ffb7
BaryFVMNorms1=fourconv(BaryFVMTest1;maxlev)

# ╔═╡ 254b0749-55bb-49fb-8a3e-de4ace4560f0
function BaryFVMTest2(grid,Λ)
# Λ is part of edge factors e!

function diffusion2(f,u,uold,sys,Δt)
	e=sys.e
	ω=sys.ω
	tris::Matrix{Int64}=sys.grid[CellNodes]
	coord::Matrix{Float64}=sys.grid[Coordinates]
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
			f[k]+=ω[inode,itri]*∇Λ∇(finitebell,coord[:,k],Λ)
		end
	end
 homogeneous_dirichlet(f,u,sys.grid)	
	
end	
    sys=EvolutionSystem(grid,diffusion2;jac=stdsparse(grid),Λ=Λ);
	
	statsolve(sys,map(finitebellx,grid))
	
end

	

# ╔═╡ 067d2d91-8f8f-4e67-90cd-e3f45c082238
fourtests(BaryFVMTest2,"baryfvm2.png")

# ╔═╡ 7545922b-6af4-4a76-bfca-10896fbf97d2
BaryFVMNorms2=fourconv(BaryFVMTest2;maxlev)

# ╔═╡ daf2aeea-e646-4364-9876-63a96b993220
md"""
### Nonlinear BaryFVM, better upwind
"""

# ╔═╡ 60e89983-7236-4653-b3ae-b400034921f5
function BaryFVMTest3(grid,Λ)


function diffusion3(f,u,uold,sys,Δt)
	e=sys.e
	ω=sys.ω
	tris::Matrix{Int64}=sys.grid[CellNodes]
	coord::Matrix{Float64}=sys.grid[Coordinates]
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
			f[k]+=ω[inode,itri]*∇Λ∇(finitebell,coord[:,k],Λ)
		end
	end
 homogeneous_dirichlet(f,u,sys.grid)	
end 
	
	sys=EvolutionSystem(grid,diffusion3;jac=stdsparse(grid),Λ=Λ);
	
	statsolve(sys,map(finitebellx,grid))
	
end



	

# ╔═╡ 2cb683c4-6050-4bc7-9ee1-e90b7b0b8117
fourtests(BaryFVMTest3,"baryfvm3.png")

# ╔═╡ 519a182b-a1ab-4e02-929e-56f2a5e3f7c4
BaryFVMNorms3=fourconv(BaryFVMTest3;maxlev)

# ╔═╡ 077b2362-7963-462e-8509-fe1da48d6cdb
md"""
### Nonlinear BaryFVM, Quenjel scheme
"""

# ╔═╡ 32e2907d-ea3a-48cc-8931-901ad15b3c57
md"""
Cut square root function
"""

# ╔═╡ 1a6b3bca-ebeb-49d8-b7af-04fc5394536f
xsqrt(x)= x<0 ? 0 : sqrt(x+1.0e-20)

# ╔═╡ 0ebac4aa-eb69-4cd3-a2eb-e6276ba07ded
function BaryFVMTest4(grid,Λ)

function diffusion4(f,u,uold,sys,Δt)
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
			f[k]+=ω[inode]*∇Λ∇(finitebell,coord[:,k],Λ)
		end
	end
 homogeneous_dirichlet(f,u,sys.grid)	
end	
	
	sys=EvolutionSystem(grid,diffusion4;jac=stdsparse(grid),Λ=Λ);
	
	statsolve(sys,map(finitebellx,grid))
	
end



	

# ╔═╡ e73007ea-ea80-456d-afc0-a28c5f6d8f06
fourtests(BaryFVMTest4,"baryfvm4.png")

# ╔═╡ 8c48a02f-fdc0-473e-8cc2-bf47fdd120c9
BaryFVMNorms4=fourconv(BaryFVMTest4;maxlev)

# ╔═╡ c651e025-eb3f-4e80-9aac-2265ffa52439
l2(case,norms)=[norms[case][i][1] for i=1:length(norms[case])]

# ╔═╡ fac87491-af52-4fbd-98b4-b0480a26da20
h1(case,norms)=[norms[case][i][2] for i=1:length(norms[case])]

# ╔═╡ 18c38fb8-e285-4ed8-a7d9-144cb4e9ab53
function plotnorms(nrm; addplot=()->())
	clf()
	function plots(case)
	xlabel("h")
	ylabel("norms")
    PyPlot.grid()
	loglog(H,nrm(case,VoronoiFVMNorms),"+-",label="Voronoi",markersize=8)
	loglog(H,nrm(case,FEMNorms),"o-",label="FEM",markersize=5)
	loglog(H,nrm(case,BaryFVMNorms1),label="Bary1")
	loglog(H,nrm(case,BaryFVMNorms2),"o-",label="Bary2")
	loglog(H,nrm(case,BaryFVMNorms3),label="Bary3")
	loglog(H,nrm(case,BaryFVMNorms4),label="Bary4")
    addplot()
    legend(loc="lower right")
	end

	subplot(221)
	title("rgrid,ΛMatrix(100,0)")
	plots(1)

	subplot(223)
	title("rgrid,ΛMatrix(100,π/4)")
	plots(2)

	subplot(222)
	title("tgrid,ΛMatrix(100,0)")
	plots(3)

	subplot(224)
	title("tgrid,ΛMatrix(100,π/4)")
	plots(4)

	
	fig=gcf()
	tight_layout()
	fig.set_size_inches(7,7)
	fig
end

# ╔═╡ 8d520cbd-08e3-41e9-93fd-4187e7b91153
let
	function addplot()
     	loglog(H,H.^2/3,color=(0,0,0),label="O(h^2)")		
    	loglog(H,H,"--",color=(0,0,0),label="O(h)")		
	end
	plotnorms(l2;addplot)
	mysave("l2conv.png")
	gcf()
end

# ╔═╡ c823d4ac-9d79-40ba-8203-6a96de0727cc
let
	function addplot()
     	loglog(H,H.^2/3,color=(0,0,0),label="O(h^2)")		
    	loglog(H,2H,"--",color=(0,0,0),label="O(h)")		
	end
	plotnorms(h1;addplot)
    mysave("h1conv.png")	
	gcf()
end

# ╔═╡ 23273fba-31cf-4f3f-b70a-866dac8d8085
html"<hr>"

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═de8b10df-42d6-43a5-b9d7-59bf7a4d99f3
# ╠═223d2526-9e76-43c4-8ba6-78b56342b56c
# ╟─9eb0b277-5a96-436e-a233-051edfa7cb24
# ╟─acaed9cb-5ed0-4137-b65c-904b26b4848f
# ╟─207a355d-ff8f-4404-90db-daf9341f4f98
# ╟─33755343-8709-4c88-95f9-9758f1d61f26
# ╟─e5a25fb2-4525-4ca4-b155-b00f90fea311
# ╟─6b417d68-7428-4ebe-8932-eec7424f100c
# ╟─cfc41a2a-3df9-4fc2-82ec-2769b25e9d56
# ╟─093cada5-9eba-40e1-bff7-d7aee5cf6b19
# ╟─d64a304c-e970-4c1a-b2d6-8314bd69df98
# ╟─9669a6e5-0593-49ed-9fa6-35a1e3986548
# ╠═344b51fb-2a23-4620-ac30-4129e7e48df7
# ╟─04838dde-4924-4b67-b9b2-1ecd534fe861
# ╠═86c83f88-bf48-48c9-ba57-b6dae3574c18
# ╠═101e6962-1046-4f98-a363-a666c026b95e
# ╠═06abe474-3862-440e-9e02-59fb897b1562
# ╟─619675e5-abea-48e1-83ef-3f0ff4395eed
# ╟─54af9427-239b-4cc2-adaf-c443edb90df3
# ╟─250cbba5-334a-494c-9d03-f72e8fb38937
# ╟─3e191aec-1121-4584-99e9-953c894fb13c
# ╟─7b273b52-be4f-4a32-9561-d4abfb79761b
# ╠═7738bc89-efd4-4b37-b5b4-0de400657f7e
# ╠═a20f6356-f05f-4927-a0dd-6226b33c1433
# ╠═ed8a0ecd-5563-4351-8291-b35cabfc1b03
# ╠═687f9d1e-4dc1-4791-a4c4-28f7e6b46983
# ╠═50809264-3603-41d9-923c-5a163f221344
# ╟─34d148a4-8994-4f7f-9774-40a63644bcee
# ╠═d86eef8a-99e1-4fc4-ba9f-213a41b86598
# ╠═4b439328-18d9-42fc-9575-85594a342b62
# ╠═77ad9955-09bd-4042-9aad-0f94fa904ff7
# ╟─8b86f35e-ff6b-482c-b150-796da7db68ee
# ╟─49b5bc2a-96df-483d-862a-e0ee22f5e0b0
# ╠═84ebe218-4dd0-440a-91ce-f37c5d06a341
# ╠═251eeb13-2730-413f-896f-5738fda9e801
# ╠═80ca622c-a1ed-4216-b25c-428e591ca32a
# ╟─97dc31e8-a4e5-4294-b214-da3e53089cc8
# ╟─648a48f0-2af4-4b9e-80f5-1926653c02ee
# ╠═823f8838-a6ee-4ccc-8f10-27510499e01b
# ╠═5f53d081-97ca-4763-9b43-5822b97dd6a2
# ╠═d5d60909-cf40-4b0f-9626-98e2b1fb853a
# ╠═e758fc09-aa0a-4ed3-852a-faa30bc9ffb7
# ╟─7511e717-3ad2-4d3b-8445-b11b5a5bb3fc
# ╟─038dd7cb-07eb-481c-9459-29afab26c9bb
# ╠═1cfc228b-fa93-49e7-98c0-0f8d8dc45f8d
# ╠═6fd997bd-1641-4d4b-a327-bea39cea4ba0
# ╠═254b0749-55bb-49fb-8a3e-de4ace4560f0
# ╠═067d2d91-8f8f-4e67-90cd-e3f45c082238
# ╠═7545922b-6af4-4a76-bfca-10896fbf97d2
# ╟─daf2aeea-e646-4364-9876-63a96b993220
# ╠═60e89983-7236-4653-b3ae-b400034921f5
# ╠═2cb683c4-6050-4bc7-9ee1-e90b7b0b8117
# ╠═519a182b-a1ab-4e02-929e-56f2a5e3f7c4
# ╟─077b2362-7963-462e-8509-fe1da48d6cdb
# ╟─32e2907d-ea3a-48cc-8931-901ad15b3c57
# ╠═1a6b3bca-ebeb-49d8-b7af-04fc5394536f
# ╠═0ebac4aa-eb69-4cd3-a2eb-e6276ba07ded
# ╠═e73007ea-ea80-456d-afc0-a28c5f6d8f06
# ╠═8c48a02f-fdc0-473e-8cc2-bf47fdd120c9
# ╠═c651e025-eb3f-4e80-9aac-2265ffa52439
# ╠═fac87491-af52-4fbd-98b4-b0480a26da20
# ╠═18c38fb8-e285-4ed8-a7d9-144cb4e9ab53
# ╠═8d520cbd-08e3-41e9-93fd-4187e7b91153
# ╠═c823d4ac-9d79-40ba-8203-6a96de0727cc
# ╟─23273fba-31cf-4f3f-b70a-866dac8d8085
