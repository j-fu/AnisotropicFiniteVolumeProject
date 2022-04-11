function coordmatrix!(C,coord, cellnodes,k)
    spacedim=size(coord,1)
    celldim=size(cellnodes,1)
    for jj=1:celldim
        C[1,jj]=1
        for ii=1:spacedim
            C[ii+1,jj]=coord[ii,cellnodes[jj,k]]
        end
    end
end

function gradient!(G,C,factdim,I)
    clu=lu!(C)
    vol=abs(det(clu))/factdim
    ldiv!(G,clu,I)
    return vol
end


function scalpro(G,dim,jl,il)
    s=0.0
    @inbounds @simd for k=1:dim
        s+=G[jl,k+1]*G[il,k+1]
    end
    return s
end

function stiffness!(S,dim,G)
    @inbounds for il=1:dim+1
        S[il,il]=scalpro(G,dim,il,il)
        for jl=il+1:dim+1
            S[il,jl]=scalpro(G,dim,jl,il)
            S[jl,il]=S[il,jl]
        end
    end
    return S
end



function stiffness!(S,dim,Λ,G)
    @inbounds for il=1:dim+1
        ΛGil=Λ*G[il,2:dim+1]
        S[il,il]=dot(ΛGil,G[il,2:dim+1])
        for jl=il+1:dim+1
            S[il,jl]=dot(ΛGil,G[jl,2:dim+1])
            S[jl,il]=S[il,jl]
        end
    end
    return S
end

function fenorms(u,coord,cellnodes)
    l2norm=0.0
    h1norm=0.0
    dim=size(coord,1)
    nnodes=dim+1
    factdim::Float64=factorial(dim)
    S=zeros(nnodes, nnodes) # local stiffness matrix
    C=zeros(nnodes,nnodes)  # local coordinate matrix
    G=zeros(nnodes, nnodes) # shape function gradients
    I=Matrix(Diagonal(ones(nnodes)))
    ncells=size(cellnodes,2)
    local_mass_matrix= [ 2.0 1.0 1.0; 1.0 2.0 1.0; 1.0  1.0  2.0 ]
    local_mass_matrix./=12.0
    for icell=1:ncells
		coordmatrix!(C,coord,cellnodes,icell)
        vol=gradient!(G,C,factdim,I)
        stiffness!(S,dim,G)
        for i  in 1:nnodes
            for j in 1:nnodes
                uij=u[cellnodes[j,icell]]*u[cellnodes[i,icell]]
                l2norm+=uij*vol*local_mass_matrix[j,i]
                h1norm+=uij*vol*S[j,i]
            end
        end
    end
    return (sqrt(l2norm),sqrt(h1norm));
end

Dirichlet()=1.0e30


function  fem_assemble!(A_h, # Global stiffness matrix
                        F_h, # Right hand side of FEM problem
                        grid, # Discretization grid  
                        Λ,
                        f::Tf,
                        β::Tβ # Boundary function
                        ) where {Tf,Tβ}
    
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    bfacenodes=grid[BFaceNodes]
    dim=size(coord,1)
    nnodes=dim+1
    factdim::Float64=factorial(dim)
    S=zeros(nnodes, nnodes) # local stiffness matrix
    C=zeros(nnodes,nnodes)  # local coordinate matrix
    BC=zeros(dim,dim)     # local boundary coordinate matrix
    G=zeros(nnodes, nnodes) # shape function gradients
    I=Matrix(Diagonal(ones(nnodes)))
    ncells=size(cellnodes,2)
    for icell=1:ncells
        coordmatrix!(C,coord,cellnodes,icell)
        vol=gradient!(G,C,factdim,I)
        stiffness!(S,dim,Λ,G)
        for il=1:nnodes
            i=cellnodes[il,icell]
            for jl=1:nnodes
                j=cellnodes[jl,icell]
                A_h[i,j]+=vol*(S[il,jl])
            end
	    x=coord[1,cellnodes[il,icell]]
	    y=coord[2,cellnodes[il,icell]]
	    F_h[i]+=f(x,y)*vol/(nnodes)
        end
    end    
    
    # Boundary part with penalty method
    nbfaces=size(bfacenodes,2)
    bfaceregions=grid[BFaceRegions]
    for ibface in 1:nbfaces
        for idim=1:dim
            i1=bfacenodes[idim,ibface];
	    x=coord[1,i1]
	    y=coord[2,i1]
            A_h[i1,i1]+=Dirichlet();
            F_h[i1]+=Dirichlet()*β(x,y)
        end
    end
end

