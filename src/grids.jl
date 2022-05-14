"""
    rgrid(;scale,h)

Create a rectangular based triangle grid in `[-scale,scale]` with stepsize h


Defined in $(joinpath("src",basename(@__FILE__)))

"""
function rgrid(;scale=1.0,h=0.1)
	DX=-scale:h:scale
	simplexgrid(DX,DX)
end

"""
    tgrid(;scale,h)

Create a triangular grid in `[-scale,scale]` with stepsize h using the Triangle mesh generator


Defined in $(joinpath("src",basename(@__FILE__)))

"""
function tgrid(;scale=1.0,h=0.1)
    # Create a SimplexGridBuilder structure which can collect
    # geometry information
    builder=SimplexGridBuilder(Generator=Triangulate)
    # Add points, record their numbers
    p1=point!(builder,-1*scale,-1*scale)
    p2=point!(builder,1*scale,-1*scale)
    p3=point!(builder,1*scale,1*scale)
    p4=point!(builder,-1*scale,1*scale)
    
    # Connect points by respective facets (segments)
    facetregion!(builder,1)
    facet!(builder,p1,p2)
    facetregion!(builder,2)
    facet!(builder,p2,p3)
    facetregion!(builder,3)
    facet!(builder,p3,p4)
    facetregion!(builder,4)
    facet!(builder,p4,p1)
    simplexgrid(builder,maxvolume=0.75*h^2)
end
