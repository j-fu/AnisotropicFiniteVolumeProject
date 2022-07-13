### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using AnisotropicFiniteVolumeProject
    using Revise
    using Polynomials
	
    using PlutoUI
    EL=PlutoUI.ExperimentalLayout
	using Tensors	
	using GridVisualize
	using Tensors
	using VoronoiFVM
	using ExtendableGrids
	import PlutoVista
	default_plotter!(PlutoVista)
	using LinearAlgebra
end

# ╔═╡ cfc41a2a-3df9-4fc2-82ec-2769b25e9d56
X=-1.1:0.01:1.1

# ╔═╡ 332aaa75-55ee-403b-96cf-93c9e502a5bb
g1=simplexgrid(X)

# ╔═╡ d64d2ece-0ef5-4922-9f7f-38b0e9e6ac4e
g2=simplexgrid(X,X)

# ╔═╡ fcfe06a0-5506-4874-96fc-00f6b14721db
fmytest=map(finitebell,g2)

# ╔═╡ 19bc326b-efb6-4546-981a-613555e14977
dfmytest=map( (x,y)->∇Λ∇(finitebell,[x,y]),g2)

# ╔═╡ 78f7e140-b5e3-4aa8-8396-23c0fcdcc23f
scalarplot(g2,fmytest,size=(200,200)),
scalarplot(g2,dfmytest,size=(200,200))


# ╔═╡ 2f24a9fc-6065-4d09-a047-79996c9156f8
fmytest1=map(finitebell,g1)

# ╔═╡ b528af01-0fdb-433a-a4b6-57d798bcf6ce
dfmytest1=map((x)->∇Λ∇(finitebell,x),g1)

# ╔═╡ 420455d7-956a-4f2d-98fc-49681d81bb79
let
	v=GridVisualizer(;size=(500,200))
	scalarplot!(v,g1,fmytest1)
	scalarplot!(v,g1,dfmytest1,clear=false)
	reveal(v)
end

# ╔═╡ e3d98f8e-05d2-4098-b02d-60ac69e50b42
xx=0:0.01:1

# ╔═╡ b2e197c1-616f-465e-b6cf-23b033947249
dp(x)=Tensors.gradient(finitebell_core,x)

# ╔═╡ 88b409a6-2db2-494b-aac0-741617a51e38
d2p(x)=Tensors.gradient(dp,x)

# ╔═╡ f5fa1f92-9cb2-41e2-9c6f-9c30cc5e6f9f
let
    v=GridVisualizer(size=(500,200),legend=:lb)
    scalarplot!(v,xx,finitebell_core.(xx),color=:red,label="p")
    scalarplot!(v,xx,dp.(xx),clear=false,color=:green,label="dp")
    scalarplot!(v,xx,d2p.(xx),clear=false,color=:blue,label="d2p")
	reveal(v)
end

# ╔═╡ d049319d-9a8e-4293-af5c-2671a5b95b64
laplaceflux(y,u,edge)=y[1]=u[1,1]-u[1,2]

# ╔═╡ 52ee2b7b-3115-406b-9dbc-27702edf2b0e
bcondition(y,u,node)=boundary_dirichlet!(y,u,node,value=0)

# ╔═╡ 1048e68d-5a94-463c-bf6a-e5d88df74373
sys1=VoronoiFVM.System(g1;flux=laplaceflux,bcondition,source=(y,node)->y[1]=-∇Λ∇( finitebell,node[1]),species=[1])

# ╔═╡ 69994ef7-874d-4a97-9109-19ea4b6c00ee
sol1=solve(sys1)

# ╔═╡ 846ccb48-2502-471a-941e-4cc363c784d1
scalarplot(g1,sol1[1,:],size=(500,200))

# ╔═╡ a4b52c86-73d2-482c-a790-8349de1ec5fe
sys2=VoronoiFVM.System(g2;flux=laplaceflux,bcondition,source=(y,node)->y[1]=-∇Λ∇( finitebell,[node[1],node[2]]),species=[1])

# ╔═╡ 0cef4a25-19b7-4790-adca-22daf45f31bd
sol2=solve(sys2)

# ╔═╡ 4c74e06a-2e97-4fdd-afca-6925582b1f25
scalarplot(g2,sol2[1,:],size=(300,300))

# ╔═╡ b4d3e741-dff0-4165-bd51-09fa9bc20aff
α=π/4;

# ╔═╡ 35a97451-04c8-4759-b60d-3fcb77de3d45
Λ11=1000;

# ╔═╡ dbdea373-a6af-4a4b-a682-db79589778ec
D=ΛMatrix(Λ11,α)

# ╔═╡ 9816bb5b-8619-4df4-a288-eaf467073ca5
function dflux(y,u,edge)
	normal=edge[:,1]-edge[:,2]
	normal=normal/norm(normal)
	d=norm(D*normal)
	y[1]=d*(u[1,1]-u[1,2])
end

# ╔═╡ 055bd291-f390-464e-910c-f3fe2196cb7e
DX=range(-1.1,1.1,length=100)

# ╔═╡ 12b79e34-23a4-41b5-8907-0f7e69687037
dgrid2=simplexgrid(DX,DX)

# ╔═╡ d5b8095a-197b-46c7-a97d-6bc60c9571a7
gridplot(dgrid2)

# ╔═╡ b2a0a944-5ec9-4919-aecc-9c85f14183ae
dsys2=VoronoiFVM.System(dgrid2;flux=dflux,bcondition,source=(y,node)->y[1]=-∇Λ∇( finitebell,[node[1],node[2]],D),species=[1])

# ╔═╡ f40084fa-fb00-41d7-83b9-d255c0095aad
dsol2=solve(dsys2)[1,:]

# ╔═╡ f5f42b6e-9129-48cd-930a-b0e53adb4b14
esol2=map((x,y)->finitebell(x,y),dgrid2)

# ╔═╡ 7da8c293-969d-440f-a498-8e8fa134b36c
femsol=fem_solve(dgrid2,D,(x,y)-> -∇Λ∇( finitebell,[x,y],D),(x,y)->0)

# ╔═╡ df92cc36-115a-4d9a-a5d3-5cf5f84c8bea
afvmsol=afvm_solve(dgrid2,D,(x,y)-> -∇Λ∇( finitebell,[x,y],D),(x,y)->0)

# ╔═╡ 56186905-a079-4de5-8d9a-13fa6fd87579
EL.grid([
	scalarplot(dgrid2,dsol2,size=(300,300)) scalarplot(dgrid2,femsol,size=(300,300))
	scalarplot(dgrid2,afvmsol,size=(300,300))  md""""""])

# ╔═╡ 97dc31e8-a4e5-4294-b214-da3e53089cc8
md"""
## AD assembly
"""

# ╔═╡ e256d9f9-8c60-4c53-8816-7b352b8b1aa6
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
			f[k]+=ω[inode,itri]*∇Λ∇(finitebell,coord[:,k],D)
		end
	end

	bfacenodes=sys.grid[BFaceNodes]
    nbfaces=size(bfacenodes,2)
    bfaceregions=sys.grid[BFaceRegions]
    for ibface in 1:nbfaces
        for idim=1:2
            i1=bfacenodes[idim,ibface];
	        x=coord[1,i1]
	        y=coord[2,i1]
            f[i1]+=1.0e30*u[i1]
        end
    end
end

# ╔═╡ eade662b-dfba-49ea-98fc-7e40a355ce67
adsys1=EvolutionSystem(dgrid2,diffusion1;jac=stdsparse(dgrid2),Λ=D);

# ╔═╡ 5f53d081-97ca-4763-9b43-5822b97dd6a2
finitebellx(x,y)=finitebell(x,y)+0.01*rand()

# ╔═╡ a4680d1a-ffb7-4f9e-81a2-f966eb70ff72
adsol1=statsolve(adsys1,map(finitebellx,dgrid2))

# ╔═╡ 021a1379-445a-4ca5-9546-d2b64b2f487b
scalarplot(dgrid2,adsol1,size=(300,300)) 

# ╔═╡ 1cfc228b-fa93-49e7-98c0-0f8d8dc45f8d
xlog(u)= u<1.0e-20 ? log(1.0e-20) : log(u)

# ╔═╡ 254b0749-55bb-49fb-8a3e-de4ace4560f0
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
			f[k]+=ω[inode,itri]*∇Λ∇(finitebell,coord[:,k],D)
		end
	end

	bfacenodes=sys.grid[BFaceNodes]
    nbfaces=size(bfacenodes,2)
    bfaceregions=sys.grid[BFaceRegions]
    for ibface in 1:nbfaces
        for idim=1:2
            i1=bfacenodes[idim,ibface];
	        x=coord[1,i1]
	        y=coord[2,i1]
            f[i1]+=1.0e30*u[i1]
        end
    end
end

# ╔═╡ 9ba2d677-36ff-4ff9-ba05-f55ad38e2cd0
adsys2=EvolutionSystem(dgrid2,diffusion2;jac=stdsparse(dgrid2),Λ=D);

# ╔═╡ f9fc167e-247d-4ac9-974b-6605c8ef76e3
adsol2=statsolve(adsys2,map(finitebellx,dgrid2))

# ╔═╡ 791e426f-275e-4765-b990-e40491e21680
scalarplot(dgrid2,adsol2,size=(300,300)) 

# ╔═╡ 60e89983-7236-4653-b3ae-b400034921f5
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
			f[k]+=ω[inode,itri]*∇Λ∇(finitebell,coord[:,k],D)
		end
	end

	bfacenodes=sys.grid[BFaceNodes]
    nbfaces=size(bfacenodes,2)
    bfaceregions=sys.grid[BFaceRegions]
    for ibface in 1:nbfaces
        for idim=1:2
            i1=bfacenodes[idim,ibface];
	        x=coord[1,i1]
	        y=coord[2,i1]
            f[i1]+=1.0e30*u[i1]
        end
    end
end

# ╔═╡ 5fe7531b-bdc4-4fee-af11-6845b74fdfd9
adsys3=EvolutionSystem(dgrid2,diffusion3;jac=stdsparse(dgrid2),Λ=D);

# ╔═╡ a862c1cc-9473-4437-a368-5fe0d7cb41b3
adsol3=statsolve(adsys3,map(finitebellx,dgrid2))

# ╔═╡ a3af6c17-339a-4ef0-bbe8-bee338c38a6a
scalarplot(dgrid2,adsol3,size=(300,300)) 

# ╔═╡ 1a6b3bca-ebeb-49d8-b7af-04fc5394536f
xsqrt(x)= x<0 ? 0 : sqrt(x+1.0e-20)

# ╔═╡ 0ebac4aa-eb69-4cd3-a2eb-e6276ba07ded
function diffusion4(f,u,uold,sys,Δt)
	e=sys.e
	ω=sys.ω
	tris::Matrix{Int64}=sys.grid[CellNodes]
	coord::Matrix{Float64}=sys.grid[Coordinates]
	ntri=size(e,2)
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
		ω,e=baryfactors(itri,ϕij,D,coord,tris)
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
			f[k]+=ω[inode]*∇Λ∇(finitebell,coord[:,k],D)
		end
	end

	bfacenodes=sys.grid[BFaceNodes]
    nbfaces=size(bfacenodes,2)
    bfaceregions=sys.grid[BFaceRegions]
    for ibface in 1:nbfaces
        for idim=1:2
            i1=bfacenodes[idim,ibface];
	        x=coord[1,i1]
	        y=coord[2,i1]
            f[i1]+=1.0e30*u[i1]
        end
    end
end

# ╔═╡ 59a477b8-b487-40ab-8dde-b7c26dd4857b
adsys4=EvolutionSystem(dgrid2,diffusion4;jac=stdsparse(dgrid2),Λ=D);

# ╔═╡ 5e44ba9c-0174-4a18-ad03-3549666ec898
adsol4=statsolve(adsys4,map(finitebellx,dgrid2))

# ╔═╡ e29c74af-0283-42d9-90e0-bbe974752091
scalarplot(dgrid2,adsol4,size=(300,300)) 

# ╔═╡ b87af7fb-2040-44dd-a22e-bdfc9aee388c
xnorm(u)=fenorms(u,dgrid2[Coordinates],dgrid2[CellNodes])[1]

# ╔═╡ fa5f82ed-1606-41ee-b070-65f4514a8404
xnorm(esol2-dsol2),xnorm(esol2-femsol),xnorm(esol2-afvmsol)

# ╔═╡ b9991c77-8d1f-411d-a6c3-3387be874624
xnorm(esol2-adsol1),xnorm(esol2-adsol2),xnorm(esol2-adsol3),xnorm(esol2-adsol4)

# ╔═╡ 23273fba-31cf-4f3f-b70a-866dac8d8085
html"<hr>"

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═cfc41a2a-3df9-4fc2-82ec-2769b25e9d56
# ╠═332aaa75-55ee-403b-96cf-93c9e502a5bb
# ╠═d64d2ece-0ef5-4922-9f7f-38b0e9e6ac4e
# ╠═fcfe06a0-5506-4874-96fc-00f6b14721db
# ╠═19bc326b-efb6-4546-981a-613555e14977
# ╠═78f7e140-b5e3-4aa8-8396-23c0fcdcc23f
# ╠═2f24a9fc-6065-4d09-a047-79996c9156f8
# ╠═b528af01-0fdb-433a-a4b6-57d798bcf6ce
# ╠═420455d7-956a-4f2d-98fc-49681d81bb79
# ╠═e3d98f8e-05d2-4098-b02d-60ac69e50b42
# ╠═b2e197c1-616f-465e-b6cf-23b033947249
# ╠═88b409a6-2db2-494b-aac0-741617a51e38
# ╠═f5fa1f92-9cb2-41e2-9c6f-9c30cc5e6f9f
# ╠═d049319d-9a8e-4293-af5c-2671a5b95b64
# ╠═52ee2b7b-3115-406b-9dbc-27702edf2b0e
# ╠═1048e68d-5a94-463c-bf6a-e5d88df74373
# ╠═69994ef7-874d-4a97-9109-19ea4b6c00ee
# ╠═846ccb48-2502-471a-941e-4cc363c784d1
# ╠═a4b52c86-73d2-482c-a790-8349de1ec5fe
# ╠═0cef4a25-19b7-4790-adca-22daf45f31bd
# ╠═4c74e06a-2e97-4fdd-afca-6925582b1f25
# ╠═b4d3e741-dff0-4165-bd51-09fa9bc20aff
# ╠═35a97451-04c8-4759-b60d-3fcb77de3d45
# ╠═dbdea373-a6af-4a4b-a682-db79589778ec
# ╠═9816bb5b-8619-4df4-a288-eaf467073ca5
# ╠═055bd291-f390-464e-910c-f3fe2196cb7e
# ╠═12b79e34-23a4-41b5-8907-0f7e69687037
# ╠═d5b8095a-197b-46c7-a97d-6bc60c9571a7
# ╠═b2a0a944-5ec9-4919-aecc-9c85f14183ae
# ╠═f40084fa-fb00-41d7-83b9-d255c0095aad
# ╠═f5f42b6e-9129-48cd-930a-b0e53adb4b14
# ╠═56186905-a079-4de5-8d9a-13fa6fd87579
# ╠═fa5f82ed-1606-41ee-b070-65f4514a8404
# ╠═7da8c293-969d-440f-a498-8e8fa134b36c
# ╠═df92cc36-115a-4d9a-a5d3-5cf5f84c8bea
# ╟─97dc31e8-a4e5-4294-b214-da3e53089cc8
# ╠═e256d9f9-8c60-4c53-8816-7b352b8b1aa6
# ╠═eade662b-dfba-49ea-98fc-7e40a355ce67
# ╠═5f53d081-97ca-4763-9b43-5822b97dd6a2
# ╠═a4680d1a-ffb7-4f9e-81a2-f966eb70ff72
# ╠═021a1379-445a-4ca5-9546-d2b64b2f487b
# ╠═1cfc228b-fa93-49e7-98c0-0f8d8dc45f8d
# ╠═254b0749-55bb-49fb-8a3e-de4ace4560f0
# ╠═9ba2d677-36ff-4ff9-ba05-f55ad38e2cd0
# ╠═f9fc167e-247d-4ac9-974b-6605c8ef76e3
# ╠═791e426f-275e-4765-b990-e40491e21680
# ╠═60e89983-7236-4653-b3ae-b400034921f5
# ╠═5fe7531b-bdc4-4fee-af11-6845b74fdfd9
# ╠═a862c1cc-9473-4437-a368-5fe0d7cb41b3
# ╠═a3af6c17-339a-4ef0-bbe8-bee338c38a6a
# ╠═1a6b3bca-ebeb-49d8-b7af-04fc5394536f
# ╠═0ebac4aa-eb69-4cd3-a2eb-e6276ba07ded
# ╠═59a477b8-b487-40ab-8dde-b7c26dd4857b
# ╠═5e44ba9c-0174-4a18-ad03-3549666ec898
# ╠═e29c74af-0283-42d9-90e0-bbe974752091
# ╠═b87af7fb-2040-44dd-a22e-bdfc9aee388c
# ╠═b9991c77-8d1f-411d-a6c3-3387be874624
# ╟─23273fba-31cf-4f3f-b70a-866dac8d8085
