"""
    finitebell_core(x)

Polynomial $(finitebell_core).
This polynomial has the following property:
- `p(0)=1`
- `p'(0)=1`
- `p''(0)=1`
- `p(1)=0`
- `p'(1)=0`
- `p''(1)=0`

Defined in $(joinpath("src",basename(@__FILE__)))
"""
const finitebell_core=Polynomial([1,0,0,-10,15,-6])


"""
    finitebell(x)
    finitebell(x,y)
    finitebell([x,y])

Twice differentiable function of one or two variables with finite support 
such that for `r=x` or `r=sqrt(x^2+x^2)`, `finitebell(r)=finitebell_core(r)` if `r<1` ,
otherwise, `finitebell(r)=0`. 

Defined in $(joinpath("src",basename(@__FILE__)))
"""
function finitebell end


function finitebell(x)
	if x<-1
		return zero(x)
	elseif x<0.0
		return finitebell_core(-x)
	elseif x<1
		return finitebell_core(x)
	else
		return zero(x)
	end
end


function finitebell(x,y)
	r=norm((x,y))
	r<1 ? finitebell_core(r) : 0.0
end

finitebell(x::AbstractVector)=finitebell(x[1],x[2])


"""
    d1finitebell(x)
Derivative of finitebell, calculated using automatic differentiation

Defined in $(joinpath("src",basename(@__FILE__)))
"""
d1finitebell(x::Number)=Tensors.gradient(finitebell,x)


"""
    d2finitebell(x)
Second derivative of finitebell, calculated using automatic differentiation

Defined in $(joinpath("src",basename(@__FILE__)))
"""
d2finitebell(x::Number)=Tensors.gradient(d1finitebell,x)




"""
     ∇Λ∇(u,x,Λ=I)

For a ``1\\times1`` or ``2\\times 2`` matrix, with the help of automatic
differentiation apply differential operator ``\\nabla\\cdot \\Lambda  \\nabla`` to a function
of one or two unknowns


Defined in $(joinpath("src",basename(@__FILE__)))
"""
function  ∇Λ∇ end

∇Λ∇(u,x,Λ=I)=Tensors.divergence(x->Vec((Λ*Tensors.gradient(u,x))...),Vec(x...))

∇Λ∇(u,x::Number,Λ=1)=∇Λ∇(x->u(x[1]),Vec(x),Λ*I)[1]


"""
    rotator(α)

Generate rotation matrix ``\\begin{pmatrix} \\cos(α)& -\\sin(α) \\\\ \\sin(α)& \\cos(α) \\end{pmatrix}``

Defined in $(joinpath("src",basename(@__FILE__)))

"""
rotator(α)=[cos(α) -sin(α); sin(α) cos(α)]


"""
    ΛMatrix(Λ11,α)


Generate the anisotropy matrix  ``\\begin{pmatrix} \\Lambda_{11}& 0 \\\\ 0 &1\\end{pmatrix}``
rotated by α.

Defined in $(joinpath("src",basename(@__FILE__)))

"""
function ΛMatrix(Λ11,α)
    r=rotator(α)
    r*[Λ11 0 ; 0 1]*r'
end

