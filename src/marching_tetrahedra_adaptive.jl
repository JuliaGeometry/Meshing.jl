
function vertices_mt(h::HyperRectangle, ::Type{SV}) where SV
    o = SV(h.origin...)
    w = SV(h.widths...)
    @inbounds (o,
     o.+SV(0,w[2],0),
     o.+SV(w[1],w[2],0),
     o.+SV(w[1],0,0),
     o.+SV(0,0,w[3]),
     o.+SV(0,w[2],w[3]),
     o.+w,
     o.+SV(w[1],0,w[3]))
end

@inline function vertPos(e, width, origin, vals::V, iso, eps, ::Type{VertType}) where {V, VertType}
    T = eltype(vals)

    ixs     = voxEdgeCrnrs[e]
    srcVal  = vals[ixs[1]]
    tgtVal  = vals[ixs[2]]
    a       = min(max((iso-srcVal)/(tgtVal-srcVal), eps), one(T)-eps)
    b       = one(T)-a
    c1 = voxCrnrPos(VertType)[ixs[1]]
    c2 = voxCrnrPos(VertType)[ixs[2]]

    (c1 .* b + c2.* a) .* width .+ origin
end


@inline function getVertId(e, width, vals, iso::Real, origin, vtsAry::Vector, eps::Real)

    VertType = eltype(vtsAry)

    # calculate vert position
    v = vertPos(e, width, origin, vals, iso, eps, VertType)
    push!(vtsAry, v)

    return length(vtsAry)
end

"""
    procVox(vals, iso::Real, x, y, z, nx, ny,
                    vts::Dict, vtsAry::Vector, fcs::Vector,
                    eps::Real)

Processes a voxel, adding any new vertices and faces to the given
containers as necessary.
"""
function procVox(vals, iso::Real, width, origin, vtsAry::Vector, fcs::Vector,
                 eps::Real, cubeindex)
    VertType = eltype(vtsAry)
    FaceType = eltype(fcs)
    # check each sub-tetrahedron in the voxel
    @inbounds for i = 1:6
        tIx = tetIx(i, cubeindex)
        (tIx == 0x00 || tIx == 0x0f) && continue

        e = tetTri[tIx]

        # add the face to the list
        push!(fcs, FaceType(
                    getVertId(voxEdgeId(i, e[1]), width, vals, iso, origin, vtsAry, eps),
                    getVertId(voxEdgeId(i, e[2]), width, vals, iso, origin, vtsAry, eps),
                    getVertId(voxEdgeId(i, e[3]), width, vals, iso, origin, vtsAry, eps)))

        # bail if there are no more faces
        iszero(e[4]) && continue
        push!(fcs, FaceType(
                    getVertId(voxEdgeId(i, e[4]), width, vals, iso, origin, vtsAry, eps),
                    getVertId(voxEdgeId(i, e[5]), width, vals, iso, origin, vtsAry, eps),
                    getVertId(voxEdgeId(i, e[6]), width, vals, iso, origin, vtsAry, eps)))
    end
end

function isosurface(f::Function, method::AdaptiveMarchingTetrahedra, ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{3, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2)) where {VertType, FaceType}

    ET = eltype(VertType)

    # arrays for vertices and faces
    vts = VertType[]
    fcs = FaceType[]

    # refinement queue
    refinement_queue = HyperRectangle{3,ET}[HyperRectangle{3,ET}(origin,widths)]

    val_store = Dict{NTuple{3,ET},ET}();

    @inbounds while true

        if isempty(refinement_queue)
            break
        end

        cell = pop!(refinement_queue)
        points = vertices_mt(cell, VertType)

        iso_vals = get_iso_vals(f,val_store,points)

        #Determine the index into the edge table which
        #tells us which vertices are inside of the surface
        cubeindex = method.insidepositive ? _get_cubeindex_pos(iso_vals, method.iso) : _get_cubeindex(iso_vals, method.iso)
        #interpindex = method.insidepositive ? _get_interpindex_pos(iso_vals, method.iso) : _get_interpindex(iso_vals, method.iso)

        value_interp = sum(iso_vals)*0.125
        value_true = get_iso_vals(f, val_store, center(cell))

        if (cubeindex == 0xff && value_true < 0) || (iszero(cubeindex) && value_true > 0)
            continue
        elseif minimum(cell.widths) > method.atol && !isapprox(value_interp, value_true, rtol=method.rtol, atol=method.atol)
            append!(refinement_queue, octsplit(cell))
        else
            # Find the vertices where the surface intersects the cube
            # The underlying space is non-linear so there will be error otherwise
            procVox(iso_vals, method.iso, cell.widths, cell.origin, vts, fcs, method.eps, cubeindex)
        end
    end
    vts,fcs
end
