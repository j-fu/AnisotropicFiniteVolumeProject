

∇Λ∇(u,x,D=I)=Tensors.divergence(x->Vec((D*Tensors.gradient(u,x))...),Vec(x...))

∇Λ∇(u,x::Number,D=1)=∇Λ∇(x->u(x[1]),Vec(x),D*I)[1]



const finitebell_core=Polynomial([1,0,0,-10,15,-6])

function finitebell(x)
	if x<-1
		return 0.0
	elseif x<0.0
		return finitebell_core(-x)
	elseif x<1
		return finitebell_core(x)
	else
		return 0
	end
end


function finitebell(x,y)
	r=norm((x,y))
	r<1 ? finitebell_core(r) : 0.0
end

finitebell(x::AbstractVector)=finitebell(x[1],x[2])
