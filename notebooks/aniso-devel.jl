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

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
    using Pkg
    Pkg.activate(joinpath(@__DIR__,".."))

    ENV["LC_NUMERIC"]="C"
	using Revise
	using AnisotropicFiniteVolumeProject
    using PlutoUI, PyPlot,SimplexGridFactory,ExtendableGrids,ExtendableSparse,GridVisualize,SparseArrays, Printf, HypertextLiteral,PlutoVista,Triangulate
	using LinearAlgebra
end;

# ╔═╡ 7489da38-186d-43ce-a4f0-d6969f6b3cb0
md"""
# Developing finite volume methods for Anisotropic Problems
"""

# ╔═╡ d1edd6da-55f6-11eb-2411-8f7d467a6ce3
TableOfContents(title="Contents",indent=true,aside=true, depth=4)

# ╔═╡ adb2d93e-55fe-11eb-021e-377f2cc0804a
# We use the SimplexGridBuilder from SimplexGridFactory.jl
function describe_grid()
    # Create a SimplexGridBuilder structure which can collect
    # geometry information
    builder=SimplexGridBuilder(Generator=Triangulate)
    
    # Add points, record their numbers
    p1=point!(builder,-1,-1)
    p2=point!(builder,1,-1)
    p3=point!(builder,1,1)
    p4=point!(builder,-1,1)
    
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
    builder
end

# ╔═╡ d8af7cd2-56b2-11eb-2d75-6dd0139a9525
builder=describe_grid()

# ╔═╡ 5e7f0de3-797f-452f-bdac-f74fed0b7231
md"""
We can plot the input and the possible output of the builder.
"""

# ╔═╡ 3ed78e0e-2915-4ece-a65f-ebb60ef5cef3
md"""
 The simplexgrid method creates an object of type ExtendableGrid which is
 defined in ExtendableGrids.jl. We can overwrite the maxvolume default
 which we used in `describe_grid`.
"""

# ╔═╡ d2674a6a-5600-11eb-0f25-53e931f4b907
md"""
#### Desired number of triangles
From the desired number of triangles, we can calculate a value fo the maximum area constraint passed to the mesh generator:
Desired number of triangles: $(@bind desired_number_of_triangles Slider(10:10:10000,default=50,show_value=true))
"""

# ╔═╡ 89213fd8-56b3-11eb-2095-bf28d88d8c8b

grid=simplexgrid(builder,maxvolume=4/desired_number_of_triangles)

# ╔═╡ 6c43ea80-52ea-4db0-9b9f-19fd838f19be
gridplot(grid, Plotter=PlutoVista,resolution=(300,300))

# ╔═╡ ab4ba264-560a-11eb-18df-2d42c80ce562
f(x,y)=sinpi(x)*sinpi(y)

# ╔═╡ af4ea19a-560a-11eb-2429-4bd8a3d7516e
 β(x,y)=0

# ╔═╡ cc714e7c-5bd6-11eb-32d1-8fbc1fef150b
Λ=1

# ╔═╡ ed4e7d23-afde-4ec9-b054-7b54074acd74
solution=fem_solve(grid,1,f,β)

# ╔═╡ a7de78dc-55ff-11eb-0973-3535cdbdbca3

scalarplot(grid,solution,Plotter=PlutoVista,resolution=(300,300),isolines=11,colormap=:bwr)

# ╔═╡ 1f3cd802-391d-406e-813c-56761ad53620
md"""
Define an exact solution of the homogeneous Dirichlet boundary value problem
on ``\Omega=(-1,1)\times (-1,1)``


$\begin{align}
  -\nabla \delta \cdot \nabla u&=f & \text{in}\; \Omega\\
   u &= 0 &\text{on}\;\partial\Omega 
\end{align}$


"""

# ╔═╡ acbfb736-5f4f-4b4b-aead-5586b2609785
k=1; l=1;

# ╔═╡ 387357c3-825a-4ba1-a232-4c6f4baafd7b
fexact(x,y)=sinpi(k*x)*sinpi(l*y);

# ╔═╡ 0a08276c-433c-4804-871e-549d6a3e1f6d
md"""
The right corresponding hand side is
"""

# ╔═╡ dbb750c7-9ec3-4cda-97cb-5a13a1311a0f
frhs(x,y)=(k^2+l^2)*pi^2*fexact(x,y);

# ╔═╡ c6c0326a-5bd9-11eb-1b53-110a5f149742
md"""
Run convergence test for a number of grid refinement levels
"""

# ╔═╡ 00abe95c-5bd1-11eb-2d50-196b1fc8a690
function convergence_test(;nref0=0, nref1=1,k=1,l=1,fem=false)
    allh=[]
    alll2=[]
    allh1=[]
    
    β(x,y)=0
    Λ=1.0 
    for iref=nref0:nref1
       # define the refinement level via the maximum area constraint
        area=0.1*2.0^(-2*iref)
        h=sqrt(area)
        grid=simplexgrid(builder,maxvolume=area)
        uexact=map(fexact,grid)
        
        n=num_nodes(grid)
        rhs=zeros(n)
        matrix=ExtendableSparseMatrix(n,n)
        rhs=zeros(n)
        if fem       
           fem_assemble!(matrix,rhs,grid,Λ,frhs,β)
		else
        	fvm_assemble!(matrix,rhs,grid,Λ,frhs,β)
		end
         sol=matrix\rhs
        
        if fem       
	        (l2norm,h1norm)=fenorms(uexact-sol,grid[Coordinates],grid[CellNodes])
		else
	        (l2norm,h1norm)=fvnorms(uexact-sol,grid[Coordinates],grid[CellNodes])
		end
		
        push!(allh,h)
        push!(allh1,h1norm)
        push!(alll2,l2norm)
    end
    allh,alll2,allh1
end


# ╔═╡ 0f2c90de-5bd2-11eb-29fd-bbc142f910ad
allh,alll2,allh1=convergence_test(nref0=0,nref1=6,fem=true)

# ╔═╡ 12940e68-5bd1-11eb-3386-235dee1f3db7
let
    PyPlot.clf()
    PyPlot.grid()
    PyPlot.xlabel("h")
    PyPlot.ylabel("error")
    PyPlot.loglog(allh, alll2, label="l2",color=:green,marker="o")
    PyPlot.loglog(allh, allh.^2, label="\$O(h^2)\$",color=:green)
    PyPlot.loglog(allh, allh1, label="h1",color=:red,marker="o")
    PyPlot.loglog(allh, allh, label="O(h)",color=:red)
    PyPlot.legend()
    gcf().set_size_inches(5,3)
    gcf()
end

# ╔═╡ ba03e561-c3ab-4a31-9719-fc50431b51fe
html"""<hr>"""

# ╔═╡ ecc60f86-0619-4d85-b67c-dd83db65f358
begin
    highlight(mdstring,color)= htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""
    
    macro important_str(s)    :(highlight(Markdown.parse($s),"#ffcccc")) end
    macro definition_str(s)    :(highlight(Markdown.parse($s),"#ccccff")) end
    macro statement_str(s)    :(highlight(Markdown.parse($s),"#ccffcc")) end
        
        
    html"""
    <style>
     h1{background-color:#dddddd;  padding: 10px;}
     h2{background-color:#e7e7e7;  padding: 10px;}
     h3{background-color:#eeeeee;  padding: 10px;}
     h4{background-color:#f7f7f7;  padding: 10px;}
    </style>
"""
end


# ╔═╡ Cell order:
# ╟─7489da38-186d-43ce-a4f0-d6969f6b3cb0
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─d1edd6da-55f6-11eb-2411-8f7d467a6ce3
# ╠═adb2d93e-55fe-11eb-021e-377f2cc0804a
# ╠═d8af7cd2-56b2-11eb-2d75-6dd0139a9525
# ╟─5e7f0de3-797f-452f-bdac-f74fed0b7231
# ╟─3ed78e0e-2915-4ece-a65f-ebb60ef5cef3
# ╠═89213fd8-56b3-11eb-2095-bf28d88d8c8b
# ╟─d2674a6a-5600-11eb-0f25-53e931f4b907
# ╠═6c43ea80-52ea-4db0-9b9f-19fd838f19be
# ╠═a7de78dc-55ff-11eb-0973-3535cdbdbca3
# ╠═ab4ba264-560a-11eb-18df-2d42c80ce562
# ╠═af4ea19a-560a-11eb-2429-4bd8a3d7516e
# ╠═cc714e7c-5bd6-11eb-32d1-8fbc1fef150b
# ╠═ed4e7d23-afde-4ec9-b054-7b54074acd74
# ╟─1f3cd802-391d-406e-813c-56761ad53620
# ╠═acbfb736-5f4f-4b4b-aead-5586b2609785
# ╠═387357c3-825a-4ba1-a232-4c6f4baafd7b
# ╟─0a08276c-433c-4804-871e-549d6a3e1f6d
# ╠═dbb750c7-9ec3-4cda-97cb-5a13a1311a0f
# ╟─c6c0326a-5bd9-11eb-1b53-110a5f149742
# ╠═00abe95c-5bd1-11eb-2d50-196b1fc8a690
# ╠═0f2c90de-5bd2-11eb-29fd-bbc142f910ad
# ╠═12940e68-5bd1-11eb-3386-235dee1f3db7
# ╟─ba03e561-c3ab-4a31-9719-fc50431b51fe
# ╟─ecc60f86-0619-4d85-b67c-dd83db65f358
