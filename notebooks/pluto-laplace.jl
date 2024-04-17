### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	using Pkg 
	
	Pkg.activate(joinpath(@__DIR__,".."))
    using PlutoUI
	using CairoMakie
	using GridVisualize
    using VoronoiFVM
	GridVisualize.default_plotter!(CairoMakie)
end

# ╔═╡ 3c9dece1-ef77-4599-a21a-9c7ca875621f
md"""
# Laplace problem with VoronoiFVM.jl
"""

# ╔═╡ 416e143c-ca48-4fa9-9203-5e0eff028a87
md"""
Set the number of unknowns in  x resp. y direction:
"""

# ╔═╡ 1c0710ca-99c2-487b-918b-a830003e01ad
n=20

# ╔═╡ dac3346b-3de4-42c2-b378-1f07b7a6a942
md"""
Define a rectangular discretizaton grid:
"""

# ╔═╡ 29d93a59-88b0-4245-b404-f351a20c06b4
X=collect(0:1.0/n:1)

# ╔═╡ f951176d-dcc8-4a71-bc9a-e6cd82abd5ea
 grid=VoronoiFVM.Grid(X,X)

# ╔═╡ 476f51da-e4d7-48c4-98e4-ae6ac5543b05
gridplot(grid,resolution=(300,300))

# ╔═╡ 6c46a4d9-97c1-43a1-9d05-76afe2328f2d
md"""
Define a two-point flux between Voronoi cells:
"""

# ╔═╡ d69caea6-e105-427f-86e8-53e8668ceebd
function laplace_flux!(f,u,edge)
    f[1]=u[1,1]-u[1,2]
end;

# ╔═╡ 23b4d53a-dc9f-4372-ae21-896729bec691
md"""
Create a finite volume system:
"""

# ╔═╡ 4537fa51-79ce-4d53-adda-1ecdc4e3b866
begin
	ispec=1
	physics=VoronoiFVM.Physics(flux=laplace_flux!)
	sys=VoronoiFVM.System(grid,physics)
    enable_species!(sys,ispec,[1])
    boundary_dirichlet!(sys,ispec,1,0.0)
    boundary_dirichlet!(sys,ispec,3,1.0)
end;

# ╔═╡ fccbcfd7-acfc-4336-b22b-e6a342c444a3
md"""
Solve: 
"""

# ╔═╡ 35e1216c-b9a3-43ad-ab47-53384371731a
solution=solve(unknowns(sys,inival=0),sys)

# ╔═╡ 6c541326-3a86-4513-b81f-18adeeb81b58
scalarplot(grid,solution[ispec,:])

# ╔═╡ 8d35fd9b-df03-441d-91c6-aadcc2abb6f8
html"""<hr>"""

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─3c9dece1-ef77-4599-a21a-9c7ca875621f
# ╟─416e143c-ca48-4fa9-9203-5e0eff028a87
# ╠═1c0710ca-99c2-487b-918b-a830003e01ad
# ╟─dac3346b-3de4-42c2-b378-1f07b7a6a942
# ╠═29d93a59-88b0-4245-b404-f351a20c06b4
# ╠═f951176d-dcc8-4a71-bc9a-e6cd82abd5ea
# ╠═476f51da-e4d7-48c4-98e4-ae6ac5543b05
# ╟─6c46a4d9-97c1-43a1-9d05-76afe2328f2d
# ╠═d69caea6-e105-427f-86e8-53e8668ceebd
# ╟─23b4d53a-dc9f-4372-ae21-896729bec691
# ╠═4537fa51-79ce-4d53-adda-1ecdc4e3b866
# ╟─fccbcfd7-acfc-4336-b22b-e6a342c444a3
# ╠═35e1216c-b9a3-43ad-ab47-53384371731a
# ╠═6c541326-3a86-4513-b81f-18adeeb81b58
# ╟─8d35fd9b-df03-441d-91c6-aadcc2abb6f8
