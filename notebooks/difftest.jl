### A Pluto.jl notebook ###
# v0.18.0

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

# ╔═╡ e47a877e-004f-4666-9caf-b0d72e25949b
using HypertextLiteral

# ╔═╡ c36e1588-4144-41c8-befc-1708e84adc2e
Revise.revise(AnisoFV)

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
dp(x)=Tensors.gradient(AnisoFV.finitebell_core,x)

# ╔═╡ 88b409a6-2db2-494b-aac0-741617a51e38
d2p(x)=Tensors.gradient(dp,x)

# ╔═╡ f5fa1f92-9cb2-41e2-9c6f-9c30cc5e6f9f
let
    v=GridVisualizer(size=(500,200),legend=:lb)
    scalarplot!(v,xx,AnisoFV.finitebell_core.(xx),color=:red,label="p")
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

# ╔═╡ b552ae4a-75e5-464b-aa82-839134531785
rot(α)=[cos(α) -sin(α); sin(α) cos(α)]

# ╔═╡ 055bd291-f390-464e-910c-f3fe2196cb7e
DX=range(-1.1,1.1,length=100)

# ╔═╡ 12b79e34-23a4-41b5-8907-0f7e69687037
dgrid2=simplexgrid(DX,DX)

# ╔═╡ f5f42b6e-9129-48cd-930a-b0e53adb4b14
esol2=map((x,y)->finitebell(x,y),dgrid2)

# ╔═╡ 23273fba-31cf-4f3f-b70a-866dac8d8085
html"<hr>"

# ╔═╡ 1f3a80ed-1c10-4be1-8294-1e6c4b69d372
function myaside(x;top=10)
	u=rand(1:10000)
	@htl("""
		<style>
		
		
		@media (min-width: calc(700px + 30px + 300px)) {
			aside.plutoui-aside-wrapper$(u) {

	color: var(--pluto-output-color);
	position:fixed;
	right: 1rem;
	top: $(top)px;
	width: 400px;
	padding: 10px;
	border: 3px solid rgba(0, 0, 0, 0.15);
	border-radius: 10px;
	box-shadow: 0 0 11px 0px #00000010;
	/* That is, viewport minus top minus Live Docs */
	max-height: calc(100vh - 5rem - 56px);
	overflow: auto;
	z-index: 40;
	background-color: var(--main-bg-color);
	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);
	
			}
			aside.plutoui-aside-wrapper > div {
#				width: 300px;
			}
		}
		</style>
		
		<aside class="plutoui-aside-wrapper$(u)">
		<div>
		$(x)
		</div>
		</aside>
		
		""")
end


# ╔═╡ cf5fd00c-b13a-42ea-8703-feba95c1fdaf
myaside(md""" 
α=$(@bind sα PlutoUI.confirm(PlutoUI.TextField(default="0"))) 

Λ11=$(@bind sΛ11 PlutoUI.confirm(PlutoUI.TextField(default="1"))) 
""")

# ╔═╡ dbdea373-a6af-4a4b-a682-db79589778ec
begin
	α=eval(Meta.parse(sα))
	Λ11=eval(Meta.parse(sΛ11))
	D0=[Λ11 0 ; 0 1]
	D=rot(α)*D0*rot(α)'
end

# ╔═╡ 9816bb5b-8619-4df4-a288-eaf467073ca5
function dflux(y,u,edge)
	normal=edge[:,1]-edge[:,2]
	normal=normal/norm(normal)
	d=norm(D*normal)
	y[1]=d*(u[1,1]-u[1,2])
end

# ╔═╡ b2a0a944-5ec9-4919-aecc-9c85f14183ae
dsys2=VoronoiFVM.System(dgrid2;flux=dflux,bcondition,source=(y,node)->y[1]=-∇Λ∇( finitebell,[node[1],node[2]],D),species=[1])

# ╔═╡ f40084fa-fb00-41d7-83b9-d255c0095aad
dsol2=solve(dsys2)[1,:]

# ╔═╡ 7da8c293-969d-440f-a498-8e8fa134b36c
femsol=fem_solve(dgrid2,D,(x,y)-> -∇Λ∇( finitebell,[x,y],D),(x,y)->0)

# ╔═╡ df92cc36-115a-4d9a-a5d3-5cf5f84c8bea
afvmsol=afvm_solve(dgrid2,D,(x,y)-> -∇Λ∇( finitebell,[x,y],D),(x,y)->0)

# ╔═╡ 56186905-a079-4de5-8d9a-13fa6fd87579
EL.grid([
	scalarplot(dgrid2,dsol2,size=(300,300)) scalarplot(dgrid2,femsol,size=(300,300))
	scalarplot(dgrid2,afvmsol,size=(300,300))  md""""""])

# ╔═╡ fa5f82ed-1606-41ee-b070-65f4514a8404
norm(esol2-dsol2),norm(esol2-femsol),norm(esol2-afvmsol)

# ╔═╡ 6aa59640-4de2-494c-bb45-6fcd6b62fa4b
myaside(htl"""D0=$(repr(D0)) <br> Λ=$(repr(round.(D,digits=3)))""",top=100)

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═c36e1588-4144-41c8-befc-1708e84adc2e
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
# ╠═b552ae4a-75e5-464b-aa82-839134531785
# ╠═cf5fd00c-b13a-42ea-8703-feba95c1fdaf
# ╠═dbdea373-a6af-4a4b-a682-db79589778ec
# ╠═6aa59640-4de2-494c-bb45-6fcd6b62fa4b
# ╠═9816bb5b-8619-4df4-a288-eaf467073ca5
# ╠═055bd291-f390-464e-910c-f3fe2196cb7e
# ╠═12b79e34-23a4-41b5-8907-0f7e69687037
# ╠═b2a0a944-5ec9-4919-aecc-9c85f14183ae
# ╠═f40084fa-fb00-41d7-83b9-d255c0095aad
# ╠═f5f42b6e-9129-48cd-930a-b0e53adb4b14
# ╠═56186905-a079-4de5-8d9a-13fa6fd87579
# ╠═fa5f82ed-1606-41ee-b070-65f4514a8404
# ╠═7da8c293-969d-440f-a498-8e8fa134b36c
# ╠═df92cc36-115a-4d9a-a5d3-5cf5f84c8bea
# ╟─23273fba-31cf-4f3f-b70a-866dac8d8085
# ╠═e47a877e-004f-4666-9caf-b0d72e25949b
# ╠═1f3a80ed-1c10-4be1-8294-1e6c4b69d372
