### A Pluto.jl notebook ###
# v0.19.41

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

# ╔═╡ 3a6b6f6d-58a8-47b4-ae6c-fc13905274fa
begin 
    using Pkg
    Pkg.activate(joinpath(@__DIR__,".."))
    ENV["LC_NUMERIC"]="C"
	using Revise
	using AnisotropicFiniteVolumeProject
	using StaticArrays
    using PlutoUI, PyPlot,SimplexGridFactory,ExtendableGrids,ExtendableSparse,GridVisualize,SparseArrays, Printf, HypertextLiteral,PlutoVista,Triangulate
	using VoronoiFVM
	using LinearAlgebra
	GridVisualize.default_plotter!(PlutoVista);
end;


# ╔═╡ c644a241-3cb6-4ef1-a453-a07186c13896
TableOfContents()

# ╔═╡ 9fed8cf3-5fad-4bba-8f09-173769099382
md"""
## Donald boxes
"""

# ╔═╡ d9af76f6-9b71-43b4-b773-2e20ae5898d8
function triplot(X; bary=false)
	if bary
   center=[sum(X[1,:])/3,sum(X[2,:])/3]
else
   center=zeros(2)
    Triangulate.tricircumcenter!(center,X[:,1],X[:,2],X[:,3])
end
	
	clf()
	
	PyPlot.fill(X[1,:],X[2,:],fill=false)
	fig=PyPlot.gcf()
	ax=PyPlot.gca()
	ax.set_aspect(1.0)
	fig.set_size_inches(4,4)
	scatter(X[1,:],X[2,:])
	scatter(center...)
	for i=1:3
		j=i%3+1
		ecenter=[X[1,i]+X[1,j],X[2,i]+X[2,j]]/2
		scatter(ecenter...)
		PyPlot.plot([center[1],ecenter[1]],[center[2],ecenter[2]])
		l=norm(center-ecenter)
		normal0=(center-ecenter)/l
		mid=(center+ecenter)/2
		normal=[normal0[2],-normal0[1]]*0.1
		scatter(mid...)
		PyPlot.plot([mid[1],mid[1]+normal[1]],[mid[2],mid[2]+normal[2]])
	end
	fig
end

# ╔═╡ 4e481553-cded-443f-a8d3-4f4cd445ccac
md"""
### Isotropic
"""

# ╔═╡ e42699e3-7ae0-4cf4-981d-ac2811b85621
function vorofactors(itri,pts,tris)
	G,vol=femgrad(itri,pts,tris)
	ω=@MVector zeros(3)
	e=@MVector zeros(3)
	AnisoFV.trifactors!(ω,e,itri,pts,tris)
	ω.=vol/3
ω,e
end

# ╔═╡ 459fa50e-6295-4662-a735-ba9f686cb608
function vfvmfactors(itri,pts,tris)
	ω=@MVector zeros(3)
	e=@MVector zeros(3)
	VoronoiFVM.cellfactors!(Triangle2D,Cartesian2D,pts,tris,itri,ω,e)
ω,e
end

# ╔═╡ e8baf009-9546-45fb-8e19-31368a2f347d
function femfactors!(ω,e,itri, coord, pointlist)
	G,vol=femgrad(itri,coord,pointlist)
	e[1]=-dot(G[2,:],G[3,:])*vol
	e[2]=-dot(G[1,:],G[3,:])*vol
	e[3]=-dot(G[1,:],G[2,:])*vol
	ω[1]=ω[2]=ω[3]=vol/3
	ω,e
end

# ╔═╡ 6aaf5454-0cf9-444a-995e-bb7bbbb275a8
femfactors(itri, coord, pointlist)=femfactors!(MVector{3,Float64}(0,0,0),MVector{3,Float64}(0,0,0),itri,coord,pointlist)

# ╔═╡ 5507f5f8-6c2e-4316-9cf5-bebbe806d987
function genfactors!(ω,e,itri,coord, pointlist;bary=true)
    G,vol=femgrad(itri,coord,pointlist)	
	
	i1=pointlist[1,itri]		
	i2=pointlist[2,itri]		
	i3=pointlist[3,itri]	


#	@info ori
if bary
   center=@SVector [sum(coord[1,pointlist[:,itri]])/3,sum(coord[2,pointlist[:,itri]])/3]

else
   center=@MVector zeros(2)
   Triangulate.tricircumcenter!(center,coord[:,i1],coord[:,i2],coord[:,i3])
   		
   	
end
 	#we assure that normal and edge vectors point in the same direction
	#Seems to work well if reference point is in the triangle, 
	#With circumcenter this does not work.
	
	ec23=@SVector [ (coord[1,i2]+coord[1,i3])/2, (coord[2,i2]+coord[2,i3])/2 ]
	e23=coord[:,i2]-coord[:,i3]
	g=(center-ec23)
	nn23=@SVector [g[2],-g[1]]
	nn23*=sign(dot(e23,nn23))
	
	ec31=@SVector [(coord[1,i1]+coord[1,i3])/2, (coord[2,i1]+coord[2,i3])/2 ]
	e31=coord[:,i3]-coord[:,i1]
	g=(center-ec31)
	nn31=@SVector [g[2],-g[1]]
	nn31*=sign(dot(e31,nn31))

	ec12=@SVector [(coord[1,i1]+coord[1,i2])/2, (coord[2,i1]+coord[2,i2])/2 ]
    e12=coord[:,i1]-coord[:,i2]
	g=(center-ec12)
	nn12=@SVector [g[2],-g[1]]
	nn12*=sign(dot(e12,nn12))

# Frolkovic's lambda values are projections of the diffusion tensor

	if false
		h23=norm(e23)
		Γ23=norm(center-ec23)
		h31=norm(e31)
		Γ31=norm(center-ec31)
		h12=norm(e12)
		Γ12=norm(center-ec12)

	
		λ12= -h12*dot(nn12-nn31,G[2,:])/Γ12
		λ23= -h23*dot(nn23-nn12,G[3,:])/Γ23
		λ31= -h31*dot(nn31-nn23,G[1,:])/Γ31
	end
# Lambda values as form factors	
	λ12= dot(nn31-nn12,G[2,:])
	λ23= dot(nn12-nn23,G[3,:])
	λ31= dot(nn23-nn31,G[1,:])
	e.=(λ23,λ31,λ12)	

	if true
    	ω[1]=ω[2]=ω[3]=vol/3
	else
	ω[1]=(h31*Γ31 + h12*Γ12)/4
	ω[2]=(h23*Γ23 + h12*Γ12)/4		
	ω[3]=(h31*Γ31 + h23*Γ23)/4
	end
	ω,e
	
end

# ╔═╡ 75a4c248-84cc-4302-bc59-51f09a6e2bb8
baryfactors(itri, coord, pointlist)=genfactors!(MVector{3,Float64}(0,0,0),MVector{3,Float64}(0,0,0),itri,coord,pointlist;bary=true)

# ╔═╡ 4fc48d46-ba39-46e6-bc14-462ec6c7ae62
circumfactors(itri, coord, pointlist)=genfactors!(MVector{3,Float64}(0,0,0),MVector{3,Float64}(0,0,0),itri,coord,pointlist;bary=false)

# ╔═╡ 55ed5e5c-b51d-4d85-afbd-677576f27f1a
begin
X=[0 1 0.0;
   0 0 1]
X=rand(2,3)
end

# ╔═╡ b185ee05-c981-4a96-a27b-849b572ce03c
femfactors(1,X,[1,2,3])

# ╔═╡ 8e05718d-33d4-4272-b9be-781ac2fa6762
baryfactors(1,X,[1,2,3])

# ╔═╡ 4fa921db-7af4-479a-aaeb-df40844f7401
circumfactors(1,X,[1,2,3])

# ╔═╡ e5523355-be26-4e47-bead-45de87a98408
vorofactors(1,X,[1,2,3])

# ╔═╡ befab178-78c3-4654-8cac-a93dbcd55c0c
vfvmfactors(1,X,[1,2,3])

# ╔═╡ baa513cd-280a-40e8-892a-7455a0804d46
triplot(X;bary=false)

# ╔═╡ 4a8b1210-69f3-49e5-bdef-d6841aafd3f3
md"""
### Anisotropic
"""

# ╔═╡ dac18f3c-9629-4f70-9680-e54cdc0d8f50
begin Y=[0 1 0.7;
   0 0 0.75]
Y=rand(2,3)
end

# ╔═╡ e259d59e-b152-4a61-af6c-e9cb2f6bc31a
Λ=[1.0 -0.1; -0.1 10.0]

# ╔═╡ 420d71cf-5aed-4f83-90c1-72dcc67ecd8b
femfactors(1,Λ,Y,[1,2,3])

# ╔═╡ c7d6215c-afb6-4513-9da3-0d9e6aa2e4c8
baryfactors(1,Λ,Y,[1,2,3])

# ╔═╡ b2581a88-e554-45e5-a691-95a14a90a541
md"""
Indeed we see that we get the same data in FEM, VFVM,DFVM
"""

# ╔═╡ e73ae183-c319-4200-b4b8-6232c28e145c
md""" 
### VFVM Test

```math
\partial_t u - \nabla \cdot (φ(u)Λ\nabla u - f(u)\mathbf V)=0
```

"""

# ╔═╡ 15cfb482-0ab3-4fd9-9395-3c4ea25c9c4e
function QSystem(nx,φ::Tφ,f::Tf,V::TV,Λ) where {Tφ,Tf,TV}
	X=0:1/nx:1
	grid=simplexgrid(X,X)
	function flux(y,u,edge)
		vh=project(edge,V)
		y[1]= φ((u[1,1]+u[1,2])/2)*Λ*(u[1,1]-u[1,2])
	    y[1]+= ( vh>0 ? f(u[1,1])*vh : f(u[1,2])*vh )
	end
    storage(y,u,node)= y[1]=u[1]
	VoronoiFVM.System(grid;flux,storage,species=[1])
end

# ╔═╡ b667fc53-cfe6-4a32-b657-23ac83f59398
sys1=QSystem(10,(u)->u^2,(u)->u^2, [1,1], 1.0e-1)

# ╔═╡ 5d274586-d7a1-46ea-b122-10f8169ff173
tsol1=solve(sys1;inival=unknowns(sys1,inival=1),times=(0,1));

# ╔═╡ cded3814-4d73-4243-91d0-8b059ae1e438
@bind t1 Slider(0:0.01:1,default=0,show_value=true)

# ╔═╡ d23a2ed0-29b8-4f7b-9e96-5c326c716de7
scalarplot(sys1.grid,tsol1(t1)[1,:],Plotter=PlutoVista)

# ╔═╡ 2a41c389-3205-4ca0-9801-7ed83bf97f89
md"""
### Generic evolution system
"""

# ╔═╡ a6c36d57-4f57-4703-9cb4-dd5a8ecd61cb
grid2=simplexgrid(-10:0.1:10)

# ╔═╡ eff86332-7395-41ab-b87f-89320511dac1
function diffusion_step(f,u,uold,grid,Δt)
	N=num_nodes(grid)
	X=view(grid[Coordinates],1,:)
	fill!(f,0.0)
	for i=1:N-1
		h=X[i+1]-X[i]
		du=(u[i]^2-u[i+1]^2)/h
		f[i]+=10*du
		f[i+1]-=10*du
		f[i]+=0.5*h*(u[i]-uold[i])/Δt
		f[i+1]+=0.5*h*(u[i+1]-uold[i+1])/Δt
	end
end

# ╔═╡ d9bd7cb2-283b-4bae-aeda-65e321bce1a4
sys2=EvolutionSystem(grid2,diffusion_step);

# ╔═╡ a5bed313-9638-47bc-99e4-2aadc521a722
tsol2=evolve(sys2,map(x->exp(-x^2),grid2));

# ╔═╡ 8b0dd553-d8e6-4310-ac45-6a4ba7c7665a
@bind t2 Slider(0:0.1:1,show_value=true)

# ╔═╡ 6004f917-7e2e-4860-a350-beb050b1241b
vis2=GridVisualizer(Plotter=PlutoVista)

# ╔═╡ abe662bd-8682-4529-8094-d1b8fb383840
scalarplot!(vis2,grid2,tsol2(t2),show=true,flimits=(0,1))

# ╔═╡ eab6d7f9-716b-4efc-a1b9-ab72c9efdc4d
md"""
### 2D Evolution system
"""

# ╔═╡ 60b6e5b2-f51c-4b81-b4c8-25b800140641
function diffusion_step2(f,u,uold,grid,Δt)
	N=num_nodes(grid)
	coord=grid[Coordinates]
	tris=grid[CellNodes]
	ntri=num_cells(grid)
	en=[2 3 ; 3 1 ; 1 2 ]'
	fill!(f,0.0)
	for itri=1:ntri
		ω,e=baryfactors(itri,coord,tris)
		for iedge=1:3
			k=tris[en[1,iedge],itri]		
			l=tris[en[2,iedge],itri]	
			du=10*e[iedge]*(u[k]-u[l])
			f[k]+=du
			f[l]-=du
		end
		for inode=1:3
			k=tris[inode,itri]
			f[k]+=ω[inode]*(u[k]-uold[k])/Δt
		end
	end
end

# ╔═╡ 64569277-bf18-48ca-831e-42c1a9f14a4a
@bind run3 CheckBox()

# ╔═╡ 8c207390-ec90-402e-9163-16c3a4660852
grid3=let
	b3=SimplexGridBuilder(Generator=Triangulate)
	p1=point!(b3,-10,-10)
	p2=point!(b3,10,-10)
	p3=point!(b3,10,10)
	p4=point!(b3,-10,10)
	facet!(b3,p1,p2)
	facet!(b3,p2,p3)
	facet!(b3,p3,p4)
	facet!(b3,p4,p1)
	g=simplexgrid(b3,maxvolume=0.1)
#	cfo=g[CellFaceOrientations]
	g
end


# ╔═╡ b64b988f-7b6e-4967-a477-778a0e72977c
gridplot(grid3,Plotter=PlutoVista);

# ╔═╡ df7c62dc-697f-41bf-90c6-fd3c8f107d1b
sys3=EvolutionSystem(grid3,diffusion_step2,jac=stdsparse(grid3));

# ╔═╡ 4ee7a84b-6d1d-49e1-bced-4811d6c3873f
inival3=map((x,y)->exp(-x^2-y^2),grid3);

# ╔═╡ 861ba01f-ea16-45b7-b00d-835950924e88
gridplot(grid3,Plotter=PlutoVista);

# ╔═╡ bd60f8d9-ef7b-4d27-ae34-e4e0eaf519d3
if run3 
	tsol3=evolve(sys3,inival3);
end

# ╔═╡ 1d9a303d-7f57-49da-82da-e061cf284c4f
@bind t3 Slider(0:0.1:1,show_value=true)

# ╔═╡ f14f58a7-8a1a-427c-a25e-23013c19109e
vis3=GridVisualizer(Plotter=PlutoVista,dim=2)

# ╔═╡ 2ad20ae2-ec6e-4df4-96de-66aeded0c8db
scalarplot!(vis3,grid3,tsol3(t3),show=true,flimits=(1,-1))

# ╔═╡ 6af270f2-9cb5-43db-8fce-9e45c6235cf3
flux3(y,u,edge)= y[1]=10*(u[1,1]-u[1,2])

# ╔═╡ 24567457-1475-4f92-aacf-46000d906883
stor3(y,u,node)= y[1]=u[1]

# ╔═╡ 20215172-ba81-44c3-afb2-2fce130ca347
vsys3=VoronoiFVM.System(grid3; flux=flux3,storage=stor3,species=[1])

# ╔═╡ fa779db9-2f16-49f2-a970-939bc47d7ac3
vtsol3=solve(vsys3,inival=reshape(inival3,1,length(inival3)),times=(0,1),
Δt=0.1,Δt_max=0.100000001,Δt_min=0.1,Δu_opt=10);

# ╔═╡ f2ab05ca-180c-45f6-a5d7-7f539af5824d
vvis3=GridVisualizer(Plotter=PlutoVista,dim=2)

# ╔═╡ f7ac95b2-b4bf-454f-8abe-fd0b0649da43
scalarplot!(vvis3,grid3,vtsol3(t3)[1,:],show=true,flimits=(1,-1))

# ╔═╡ 7318ae09-5301-429e-ad16-9a2083f815a2
norm(tsol3(t3)-vtsol3(t3)[1,:])

# ╔═╡ d35efee7-e5d3-46ab-8b94-141a5c5e8041
md"""
### 2D Anisotropic
"""

# ╔═╡ 5b627845-30e7-4502-aff8-11a88f3a0f0d
Λ4=[1 0 ; 0 1000]

# ╔═╡ 017de554-19ed-4dbd-9c76-a198207c9ec5
u4ex(x,y,t)=0.5*(cospi(x)*exp(-π^2*Λ4[1,1]*t)+1)

# ╔═╡ ceb34bc4-635a-4b20-8bc4-9ee9eb05555f
grid4=let
	b3=SimplexGridBuilder(Generator=Triangulate)
	p1=point!(b3,0,0)
	p2=point!(b3,1,0)
	p3=point!(b3,1,1)
	p4=point!(b3,0,1)
	facet!(b3,p1,p2)
	facet!(b3,p2,p3)
	facet!(b3,p3,p4)
	facet!(b3,p4,p1)
	simplexgrid(b3,maxvolume=0.001)
end


# ╔═╡ 5b1e8fc8-d2f0-45b9-a336-5637512bce47
vis4ex=GridVisualizer(dim=2,size=(300,300))

# ╔═╡ 43a46703-4fe4-478b-a202-5d528912e960
function diffusion_step4(f,u,uold,grid,Δt)
	N=num_nodes(grid)
	coord=grid[Coordinates]
	tris=grid[CellNodes]
	ntri=num_cells(grid)
	en=[2 3 ; 3 1 ; 1 2 ]'
	fill!(f,0.0)
	for itri=1:ntri
		ω,e=baryfactors(itri,Λ4,coord,tris)
		for iedge=1:3
			k=tris[en[1,iedge],itri]		
			l=tris[en[2,iedge],itri]	
			du=10*e[iedge]*(u[k]-u[l])
			f[k]+=du
			f[l]-=du
		end
		for inode=1:3
			k=tris[inode,itri]
			f[k]+=ω[inode]*(u[k]-uold[k])/Δt
		end
	end
end

# ╔═╡ b8392f4f-e517-4aa3-8e2a-8d5508304ac1
sys4=EvolutionSystem(grid4,diffusion_step4,jac=stdsparse(grid4));

# ╔═╡ 517ac05b-8161-4fa6-8c40-99b69ff64989
@bind t4 Slider(0:0.01:0.15,show_value=true)

# ╔═╡ af554f23-3d60-4ce5-b63e-4553b4dff4da
uex4t=map( (x,y)->u4ex(x,y,t4),grid4)

# ╔═╡ c7581e3b-a3da-45e8-96d3-b64d6a364218
scalarplot!(vis4ex,grid4,uex4t,show=true)

# ╔═╡ 020d5637-7277-429b-9588-7a85a64e00ec
tsol4=evolve(sys4,map( (x,y)->u4ex(x,y,0),grid4),tend=0.15,nsteps=100);

# ╔═╡ b9bc6777-9d54-40a5-838d-36809cb30e27
vis4=GridVisualizer(dim=2,size=(300,300))

# ╔═╡ 8b6d29ed-6f5b-4477-967b-e588c8833dce
scalarplot!(vis4,grid4,tsol4(t4),show=true)

# ╔═╡ 02c75f02-fc53-43dd-b620-5bd4a4e69fcf
norm(tsol4(t4)-uex4t)

# ╔═╡ Cell order:
# ╠═3a6b6f6d-58a8-47b4-ae6c-fc13905274fa
# ╠═c644a241-3cb6-4ef1-a453-a07186c13896
# ╟─9fed8cf3-5fad-4bba-8f09-173769099382
# ╠═d9af76f6-9b71-43b4-b773-2e20ae5898d8
# ╠═4e481553-cded-443f-a8d3-4f4cd445ccac
# ╠═e42699e3-7ae0-4cf4-981d-ac2811b85621
# ╠═459fa50e-6295-4662-a735-ba9f686cb608
# ╠═e8baf009-9546-45fb-8e19-31368a2f347d
# ╠═6aaf5454-0cf9-444a-995e-bb7bbbb275a8
# ╠═5507f5f8-6c2e-4316-9cf5-bebbe806d987
# ╠═75a4c248-84cc-4302-bc59-51f09a6e2bb8
# ╠═4fc48d46-ba39-46e6-bc14-462ec6c7ae62
# ╠═55ed5e5c-b51d-4d85-afbd-677576f27f1a
# ╠═b185ee05-c981-4a96-a27b-849b572ce03c
# ╠═8e05718d-33d4-4272-b9be-781ac2fa6762
# ╠═4fa921db-7af4-479a-aaeb-df40844f7401
# ╠═e5523355-be26-4e47-bead-45de87a98408
# ╠═befab178-78c3-4654-8cac-a93dbcd55c0c
# ╠═baa513cd-280a-40e8-892a-7455a0804d46
# ╟─4a8b1210-69f3-49e5-bdef-d6841aafd3f3
# ╠═dac18f3c-9629-4f70-9680-e54cdc0d8f50
# ╠═e259d59e-b152-4a61-af6c-e9cb2f6bc31a
# ╠═420d71cf-5aed-4f83-90c1-72dcc67ecd8b
# ╠═c7d6215c-afb6-4513-9da3-0d9e6aa2e4c8
# ╠═b2581a88-e554-45e5-a691-95a14a90a541
# ╠═e73ae183-c319-4200-b4b8-6232c28e145c
# ╠═15cfb482-0ab3-4fd9-9395-3c4ea25c9c4e
# ╠═b667fc53-cfe6-4a32-b657-23ac83f59398
# ╠═5d274586-d7a1-46ea-b122-10f8169ff173
# ╠═cded3814-4d73-4243-91d0-8b059ae1e438
# ╠═d23a2ed0-29b8-4f7b-9e96-5c326c716de7
# ╠═2a41c389-3205-4ca0-9801-7ed83bf97f89
# ╠═a6c36d57-4f57-4703-9cb4-dd5a8ecd61cb
# ╠═eff86332-7395-41ab-b87f-89320511dac1
# ╠═d9bd7cb2-283b-4bae-aeda-65e321bce1a4
# ╠═a5bed313-9638-47bc-99e4-2aadc521a722
# ╠═8b0dd553-d8e6-4310-ac45-6a4ba7c7665a
# ╠═6004f917-7e2e-4860-a350-beb050b1241b
# ╠═abe662bd-8682-4529-8094-d1b8fb383840
# ╠═eab6d7f9-716b-4efc-a1b9-ab72c9efdc4d
# ╠═60b6e5b2-f51c-4b81-b4c8-25b800140641
# ╠═64569277-bf18-48ca-831e-42c1a9f14a4a
# ╠═8c207390-ec90-402e-9163-16c3a4660852
# ╠═b64b988f-7b6e-4967-a477-778a0e72977c
# ╠═df7c62dc-697f-41bf-90c6-fd3c8f107d1b
# ╠═4ee7a84b-6d1d-49e1-bced-4811d6c3873f
# ╠═861ba01f-ea16-45b7-b00d-835950924e88
# ╠═bd60f8d9-ef7b-4d27-ae34-e4e0eaf519d3
# ╠═1d9a303d-7f57-49da-82da-e061cf284c4f
# ╠═f14f58a7-8a1a-427c-a25e-23013c19109e
# ╠═2ad20ae2-ec6e-4df4-96de-66aeded0c8db
# ╠═6af270f2-9cb5-43db-8fce-9e45c6235cf3
# ╠═24567457-1475-4f92-aacf-46000d906883
# ╠═20215172-ba81-44c3-afb2-2fce130ca347
# ╠═fa779db9-2f16-49f2-a970-939bc47d7ac3
# ╠═f2ab05ca-180c-45f6-a5d7-7f539af5824d
# ╠═f7ac95b2-b4bf-454f-8abe-fd0b0649da43
# ╠═7318ae09-5301-429e-ad16-9a2083f815a2
# ╠═d35efee7-e5d3-46ab-8b94-141a5c5e8041
# ╠═5b627845-30e7-4502-aff8-11a88f3a0f0d
# ╠═017de554-19ed-4dbd-9c76-a198207c9ec5
# ╠═ceb34bc4-635a-4b20-8bc4-9ee9eb05555f
# ╠═5b1e8fc8-d2f0-45b9-a336-5637512bce47
# ╠═af554f23-3d60-4ce5-b63e-4553b4dff4da
# ╠═c7581e3b-a3da-45e8-96d3-b64d6a364218
# ╠═43a46703-4fe4-478b-a202-5d528912e960
# ╠═b8392f4f-e517-4aa3-8e2a-8d5508304ac1
# ╠═517ac05b-8161-4fa6-8c40-99b69ff64989
# ╠═020d5637-7277-429b-9588-7a85a64e00ec
# ╠═b9bc6777-9d54-40a5-838d-36809cb30e27
# ╠═8b6d29ed-6f5b-4477-967b-e588c8833dce
# ╠═02c75f02-fc53-43dd-b620-5bd4a4e69fcf
