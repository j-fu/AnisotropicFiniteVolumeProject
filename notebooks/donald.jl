### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 3a6b6f6d-58a8-47b4-ae6c-fc13905274fa
begin 
    using Pkg
    Pkg.activate(joinpath(@__DIR__,".."))
    ENV["LC_NUMERIC"]="C"
	using Revise
	using AnisoFV
	using StaticArrays
    using PlutoUI, PyPlot,SimplexGridFactory,ExtendableGrids,ExtendableSparse,GridVisualize,SparseArrays, Printf, HypertextLiteral,PlutoVista,Triangulate
	using LinearAlgebra
end;


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

# ╔═╡ e42699e3-7ae0-4cf4-981d-ac2811b85621
function vorofactors(itri,pts,tris)
	ω=zeros(3)
	e=zeros(3)
	AnisoFV.trifactors!(ω,e,itri,pts,tris)
ω,e
end

# ╔═╡ 41e5fa28-ad78-467e-8f76-230f54164b2a
function femgrad(itri,coord,pointlist)
	C=@MMatrix zeros(3,3)
	AnisoFV.coordmatrix!(C,coord,pointlist,1)
	G0= @MMatrix zeros(3,3)
	I=@MMatrix zeros(3,3)
	for i=1:3 
		I[i,i]=1
	end
    vol=abs(det(C))/2
	G=view(C\I,:,2:3)
	G,vol
end

# ╔═╡ e8baf009-9546-45fb-8e19-31368a2f347d
function femfactors!(ω,e,itri, coord, pointlist)
	G,vol=femgrad(itri,coord,pointlist)
	e[1]=-dot(G[2,:],G[3,:])*vol
	e[2]=-dot(G[1,:],G[3,:])*vol
	e[3]=-dot(G[1,:],G[2,:])*vol
	ω,e
end

# ╔═╡ 6aaf5454-0cf9-444a-995e-bb7bbbb275a8
femfactors(itri, coord, pointlist)=femfactors!(MVector{3,Float64}(0,0,0),MVector{3,Float64}(0,0,0),itri,coord,pointlist)

# ╔═╡ 16117f28-9205-4380-b6e7-f9670ba7e110
function genfactors!(ω,e,itri, coord, pointlist;bary=false)
    G,vol=femgrad(itri,coord,pointlist)	

	i1=pointlist[1,itri]		
	i2=pointlist[2,itri]		
	i3=pointlist[3,itri]	

if bary
   center=@SVector [sum(coord[1,pointlist[:,itri]])/3,sum(coord[2,pointlist[:,itri]])/3]
else
   center=@MVector zeros(2)
   Triangulate.tricircumcenter!(center,coord[:,1],coord[:,2],coord[:,3])
end
	
	ec23=@SVector [ (coord[1,i2]+coord[1,i3])/2, (coord[2,i2]+coord[2,i3])/2 ]
	h23=norm(coord[:,i2]-coord[:,i3])
	Γ23=norm(center-ec23)
	g=(center-ec23)
	nn23=@SVector [g[2],-g[1]]

	ec31=@SVector [(coord[1,i1]+coord[1,i3])/2, (coord[2,i1]+coord[2,i3])/2 ]
	h31=norm(coord[:,i3]-coord[:,i1])
	Γ31=norm(center-ec31)
	g=(center-ec31)
	nn31=@SVector [g[2],-g[1]]

	ec12=@SVector [(coord[1,i1]+coord[1,i2])/2, (coord[2,i1]+coord[2,i2])/2 ]
	h12=norm(coord[:,i1]-coord[:,i2])
	Γ12=norm(center-ec12)
	g=(center-ec12)
	nn12=@SVector [g[2],-g[1]]

# Frolkovic's lambda values are projections of the diffusion tensor
	λ12= h12*dot(nn12-nn31,G[2,:])/Γ12
	λ23= h23*dot(nn23-nn12,G[3,:])/Γ23
	λ31= h31*dot(nn31-nn23,G[1,:])/Γ31
    λ23,λ31,λ12	

# Lambda values as form factors	
	λ12= dot(nn12-nn31,G[2,:])
	λ23= dot(nn23-nn12,G[3,:])
	λ31= dot(nn31-nn23,G[1,:])
	e.=(λ23,λ31,λ12)	

	ω,e
	
end

# ╔═╡ 75a4c248-84cc-4302-bc59-51f09a6e2bb8
baryfactors(itri, coord, pointlist)=genfactors!(MVector{3,Float64}(0,0,0),MVector{3,Float64}(0,0,0),itri,coord,pointlist;bary=true)

# ╔═╡ 4fc48d46-ba39-46e6-bc14-462ec6c7ae62
circumfactors(itri, coord, pointlist)=genfactors!(MVector{3,Float64}(0,0,0),MVector{3,Float64}(0,0,0),itri,coord,pointlist;bary=false)

# ╔═╡ 55ed5e5c-b51d-4d85-afbd-677576f27f1a
X=[0 1 0.6;
   0 0 0.3]

# ╔═╡ b185ee05-c981-4a96-a27b-849b572ce03c
femfactors(1,X,[1,2,3])

# ╔═╡ 8e05718d-33d4-4272-b9be-781ac2fa6762
baryfactors(1,X,[1,2,3])

# ╔═╡ 4fa921db-7af4-479a-aaeb-df40844f7401
circumfactors(1,X,[1,2,3])

# ╔═╡ e5523355-be26-4e47-bead-45de87a98408
vorofactors(1,X,[1,2,3])

# ╔═╡ baa513cd-280a-40e8-892a-7455a0804d46
triplot(X;bary=false)

# ╔═╡ Cell order:
# ╠═3a6b6f6d-58a8-47b4-ae6c-fc13905274fa
# ╟─9fed8cf3-5fad-4bba-8f09-173769099382
# ╠═d9af76f6-9b71-43b4-b773-2e20ae5898d8
# ╠═e42699e3-7ae0-4cf4-981d-ac2811b85621
# ╠═41e5fa28-ad78-467e-8f76-230f54164b2a
# ╠═e8baf009-9546-45fb-8e19-31368a2f347d
# ╠═6aaf5454-0cf9-444a-995e-bb7bbbb275a8
# ╠═16117f28-9205-4380-b6e7-f9670ba7e110
# ╠═75a4c248-84cc-4302-bc59-51f09a6e2bb8
# ╠═4fc48d46-ba39-46e6-bc14-462ec6c7ae62
# ╠═55ed5e5c-b51d-4d85-afbd-677576f27f1a
# ╠═b185ee05-c981-4a96-a27b-849b572ce03c
# ╠═8e05718d-33d4-4272-b9be-781ac2fa6762
# ╠═4fa921db-7af4-479a-aaeb-df40844f7401
# ╠═e5523355-be26-4e47-bead-45de87a98408
# ╠═baa513cd-280a-40e8-892a-7455a0804d46
