# Code derived from Robin Deits' RegionTrees and AdaptiveDistanceFields
# The RegionTrees.jl package is licensed under the MIT "Expat" License
# The AdaptiveDistanceFields.jl package is licensed under the MIT "Expat" License

# this is designed to avoid allocation of a full octtree and generate a mesh
# while sampling

function vertices(h::HyperRectangle)
    SV = SVector{3,Float64}
    o = SV(h.origin...)
    w = SV(h.widths...)
    (o,
     o.+SV(w[1],0,0),
     o.+SV(w[1],w[2],0),
     o.+SV(0,w[2],0),
     o.+SV(0,0,w[3]),
     o.+SV(w[1],0,w[3]),
     o.+w,
     o.+SV(0,w[2],w[3]))
end

center(rect::HyperRectangle) = rect.origin + 0.5 * rect.widths

function octsplit(h::HyperRectangle)
    hw = h.widths
    nw = h.widths./2
    no = h.origin
    (HyperRectangle(no, nw),
     HyperRectangle(no.+nw, nw),
     HyperRectangle(no.+Vec(nw[1],0,0), nw),
     HyperRectangle(no.+Vec(0,nw[2],0), nw),
     HyperRectangle(no.+Vec(0,0,nw[3]), nw),
     HyperRectangle(no.+Vec(nw[1],nw[2],0), nw),
     HyperRectangle(no.+Vec(nw[1],0,nw[3]), nw),
     HyperRectangle(no.+Vec(0,nw[2],nw[3]), nw),
     HyperRectangle(no.+Vec(0,0,nw[3]), nw))
end

function isosurface(f::Function, method::AdaptiveMarchingCubes, ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{3, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2)) where {VertType, FaceType}

    # arrays for vertices and faces
    vts = VertType[]
    fcs = FaceType[]

    # refinement queue
    refinement_queue = HyperRectangle{3,Float64}[HyperRectangle{3,Float64}(origin,widths)]

    while true

        if isempty(refinement_queue)
            break
        end

        cell = pop!(refinement_queue)
        points = vertices(cell)

        iso_vals = (f(points[1]),
                    f(points[2]),
                    f(points[3]),
                    f(points[4]),
                    f(points[5]),
                    f(points[6]),
                    f(points[7]),
                    f(points[8]))

        #Determine the index into the edge table which
        #tells us which vertices are inside of the surface
        cubeindex = method.insidepositive ? _get_cubeindex_pos(iso_vals, method.iso) : _get_cubeindex(iso_vals, method.iso)

        value_interp = sum(iso_vals)
        value_true = f(center(cell))
        if (cubeindex == 0xff && value_true < 0) || (iszero(cubeindex) && value_true > 0)
            continue
        elseif minimum(cell.widths) > method.atol && !isapprox(value_interp, value_true, rtol=method.rtol, atol=method.atol)
            append!(refinement_queue, octsplit(cell))
        else
            # Find the vertices where the surface intersects the cube
            # The underlying space is non-linear so there will be error otherwise
            vertlist = find_vertices_interp(points, iso_vals, cubeindex, method.iso, method.eps)

            # Create the triangle
            method.reduceverts && _mc_unique_triangles!(vts, fcs, vertlist, cubeindex, FaceType)
            !method.reduceverts && _mc_create_triangles!(vts, fcs, vertlist, cubeindex, FaceType)
        end
    end
    vts,fcs
end
