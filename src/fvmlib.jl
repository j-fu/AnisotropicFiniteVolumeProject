function trifactors!(ω, e, itri, pointlist, trianglelist)
    # Obtain the node numbers for triangle itri
    i1=trianglelist[1,itri]
    i2=trianglelist[2,itri]
    i3=trianglelist[3,itri]
    
    # Calculate triangle area: 
    #   Matrix of edge vectors
    V11= pointlist[1,i2]- pointlist[1,i1]
    V21= pointlist[2,i2]- pointlist[2,i1]
    
    V12= pointlist[1,i3]- pointlist[1,i1]
    
    V22= pointlist[2,i3]- pointlist[2,i1]
    V13= pointlist[1,i3]- pointlist[1,i2]
    V23= pointlist[2,i3]- pointlist[2,i2]
    
    #   Compute determinant 
    det=V11*V22 - V12*V21
   
    #   Area
    area=0.5*det
    area=abs(area)
    # Squares of edge lengths
    dd1=V13*V13+V23*V23 # l32
    dd2=V12*V12+V22*V22 # l31
    dd3=V11*V11+V21*V21 # l21
        
    # Contributions to e_kl=σ_kl/h_kl
    e[1]= (dd2+dd3-dd1)*0.125/area
    e[2]= (dd3+dd1-dd2)*0.125/area
    e[3]= (dd1+dd2-dd3)*0.125/area
    
    # Contributions to ω_k
    ω[1]= (e[3]*dd3+e[2]*dd2)*0.25
    ω[2]= (e[1]*dd1+e[3]*dd3)*0.25
    ω[3]= (e[2]*dd2+e[1]*dd1)*0.25
end                              



function fvm_assemble!(A_h, # System matrix
                       f_h,    # Right hand side vector
                       grid,
                       Λ,      # heat conduction coefficient 
                       f::Tf, # Source/sink function
                       β::Tβ  # boundary condition function
                       ) where{Tf,Tβ}
    penalty=1.0e30
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    bfacenodes=grid[BFaceNodes]
    
    num_nodes_per_cell=3;
    num_edges_per_cell=3;
    num_nodes_per_bface=2
    ntri=size(cellnodes,2)
    nbface=size(bfacenodes,2)
    
    # Local edge-node connectivity
    local_edgenodes=[ 2 3; 3 1; 1 2]'
    
    # Storage for form factors
    e=zeros(num_nodes_per_cell)
    ω=zeros(num_edges_per_cell)
    γ=zeros(num_nodes_per_bface)
    
    # Initialize right hand side to zero
    f_h.=0.0
    
    # Loop over all triangles
    for itri=1:ntri
        trifactors!(ω,e,itri,coord,cellnodes)
        # Assemble nodal contributions to right hand side
        for k_local=1:num_nodes_per_cell
            k_global=cellnodes[k_local,itri]
            x=coord[1,k_global]
            y=coord[2,k_global]
            f_h[k_global]+=f(x,y)*ω[k_local]
        end
        
        # Assemble edge contributions to matrix
        for iedge=1:num_edges_per_cell
            k_global=cellnodes[local_edgenodes[1,iedge],itri]
            l_global=cellnodes[local_edgenodes[2,iedge],itri]
            A_h[k_global,k_global]+=Λ*e[iedge]
            A_h[l_global,k_global]-=Λ*e[iedge]
            A_h[k_global,l_global]-=Λ*e[iedge]
            A_h[l_global,l_global]+=Λ*e[iedge]
        end
    end
    
    # Assemble boundary conditions
    
    for ibface=1:nbface
        for k_local=1:num_nodes_per_bface
            k_global=bfacenodes[k_local,ibface]
            A_h[k_global,k_global]+=penalty
            x=coord[1,k_global]
            y=coord[2,k_global]
            f_h[k_global]+=penalty*β(x,y)
        end
    end
end

function fvnorms(u,pointlist,trianglelist)
    local_edgenodes=[ 2 3; 3 1; 1 2]'
    num_nodes_per_cell=3;
    num_edges_per_cell=3;
    e=zeros(num_nodes_per_cell)
    ω=zeros(num_edges_per_cell)
    l2norm=0.0
    h1norm=0.0
    ntri=size(trianglelist,2)
    for itri=1:ntri
        trifactors!(ω,e,itri,pointlist,trianglelist)
        for k_local=1:num_nodes_per_cell
            k=trianglelist[k_local,itri]
            x=pointlist[1,k]
            y=pointlist[2,k]
            l2norm+=u[k]^2*ω[k_local]
        end
        for iedge=1:num_edges_per_cell
            k=trianglelist[local_edgenodes[1,iedge],itri]
            l=trianglelist[local_edgenodes[2,iedge],itri]
            h1norm+=(u[k]-u[l])^2*e[iedge]
        end
    end
    return (sqrt(l2norm),sqrt(h1norm));
end

