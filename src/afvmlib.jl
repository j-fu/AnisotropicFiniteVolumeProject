function femfactors!(ω,e,itri, Λ,coord, pointlist)
    G,vol=femgrad(itri,coord,pointlist)
    e[1]=-dot(G[2,:],Λ*G[3,:])*vol
    e[2]=-dot(G[1,:],Λ*G[3,:])*vol
    e[3]=-dot(G[1,:],Λ*G[2,:])*vol
    ω[1]=ω[2]=ω[3]=vol/3
    ω,e
end

femfactors(itri,Λ, coord, pointlist)=femfactors!(MVector{3,Float64}(0,0,0),MVector{3,Float64}(0,0,0),itri,Λ,coord,pointlist)

function genfactors!(ω,e,itri,Λ, coord, pointlist;bary=false)
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
    
    h23=norm(e23)
    Γ23=norm(center-ec23)
    h31=norm(e31)
    Γ31=norm(center-ec31)
    h12=norm(e12)
    Γ12=norm(center-ec12)
    
    
    β=@MMatrix zeros(3,3)
    β[1,1]=Γ12*Γ12*dot(Λ*nn12,nn12)
    β[1,2]=Γ12*Γ23*dot(Λ*nn12,nn23)
    β[1,3]=Γ12*Γ31*dot(Λ*nn12,nn31)
    β[2,1]=Γ23*Γ12*dot(Λ*nn23,nn12)
    β[2,2]=Γ23*Γ23*dot(Λ*nn23,nn23)
    β[2,3]=Γ23*Γ31*dot(Λ*nn23,nn31)
    β[3,1]=Γ31*Γ12*dot(Λ*nn31,nn12)
    β[3,2]=Γ31*Γ23*dot(Λ*nn31,nn23)
    β[3,3]=Γ31*Γ31*dot(Λ*nn31,nn31)
    
    
    # Frolkovic's lambda values are projections of the diffusion tensor
    λ12= h12*dot(Λ*(nn12-nn31),G[2,:])/Γ12
    λ23= h23*dot(Λ*(nn23-nn12),G[3,:])/Γ23
    λ31= h31*dot(Λ*(nn31-nn23),G[1,:])/Γ31
    λ23,λ31,λ12	
    
    # Lambda values as form factors	
    λ12= -dot(Λ*(nn12-nn31),G[2,:])
    λ23= -dot(Λ*(nn23-nn12),G[3,:])
    λ31= -dot(Λ*(nn31-nn23),G[1,:])
    e.=(λ23,λ31,λ12)	
    
    if bary
    	ω[1]=ω[2]=ω[3]=vol/3
    else
	ω[1]=(h31*Γ31 + h12*Γ12)/4
	ω[2]=(h23*Γ23 + h12*Γ12)/4		
	ω[3]=(h31*Γ31 + h23*Γ23)/4
    end
    ω,e,β
    
end

baryfactors(itri, Λ,coord, pointlist)=genfactors!(MVector{3,Float64}(0,0,0),MVector{3,Float64}(0,0,0),itri,Λ,coord,pointlist;bary=true)

circumfactors(itri, Λ,coord, pointlist)=genfactors!(MVector{3,Float64}(0,0,0),MVector{3,Float64}(0,0,0),itri,Λ,coord,pointlist;bary=false)



function afvm_assemble!(A_h, # System matrix
                       f_h,    # Right hand side vector
                        grid,
                        Λ,      # heat conduction coefficient 
                        f::Tf, # Source/sink function
                        β::Tβ  # boundary condition function
                        ) where{Tf,Tβ}
    N=num_nodes(grid)
    coord=grid[Coordinates]
    tris=grid[CellNodes]
    bfacenodes=grid[BFaceNodes]
    ntri=num_cells(grid)
    en=[2 3 ; 3 1 ; 1 2 ]'

    penalty=1.0e30
    num_nodes_per_cell=3;
    num_edges_per_cell=3;
    num_nodes_per_bface=2
    nbface=size(bfacenodes,2)
    
    
    # Initialize right hand side to zero
    f_h.=0.0


    for itri=1:ntri
	ω,e=baryfactors(itri,Λ,coord,tris)
	for iedge=1:3
	    k=tris[en[1,iedge],itri]		
	    l=tris[en[2,iedge],itri]
	    A_h[k,k]+=e[iedge]
	    A_h[k,l]-=e[iedge]
	    A_h[l,l]+=e[iedge]
	    A_h[l,k]-=e[iedge]
	end
	for inode=1:3
	    k=tris[inode,itri]
            x=coord[1,k]
            y=coord[2,k]
	    f_h[k]+=ω[inode]*f(x,y)
	end
    end

    # Boundary part with penalty method
    nbfaces=size(bfacenodes,2)
    bfaceregions=grid[BFaceRegions]
    for ibface in 1:nbfaces
        for idim=1:2
            i1=bfacenodes[idim,ibface];
	    x=coord[1,i1]
	    y=coord[2,i1]
            A_h[i1,i1]+=Dirichlet();
            f_h[i1]+=Dirichlet()*β(x,y)
        end
    end
end



function afvm_solve(grid,Λ,f,β)
    # Initialize sparse matrix and right hand side
    n=num_nodes(grid)
    matrix=spzeros(n,n)
    rhs=zeros(n)
    # Call the assemble function.
    afvm_assemble!(matrix,rhs,grid,Λ,f,β)
    sol=matrix\rhs
end 
