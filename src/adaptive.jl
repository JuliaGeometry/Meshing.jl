# Code derived from Robin Deits' RegionTrees and AdaptiveDistanceFields
# The RegionTrees.jl package is licensed under the MIT "Expat" License
# The AdaptiveDistanceFields.jl package is licensed under the MIT "Expat" License

# this is designed to avoid allocation of a full octtree and generate a mesh
# while sampling

struct HyperRectangle{N,T}
    origin::NTuple{N,T}
    widths::NTuple{N,T}
end


function vertices(h::HyperRectangle, ::Type{SV}) where SV
    z = zero(eltype(SV))
    @inbounds begin
        o = SV(h.origin[1],h.origin[2],h.origin[3])
        w = SV(h.widths[1],h.widths[2],h.widths[3])
        (o,
        o.+SV(w[1],z,z),
        o.+SV(w[1],w[2],z),
        o.+SV(z,w[2],z),
        o.+SV(z,z,w[3]),
        o.+SV(w[1],z,w[3]),
        o.+w,
        o.+SV(z,w[2],w[3]))
    end
end

function vertices_mt(h::HyperRectangle, ::Type{SV}) where SV
    z = zero(eltype(SV))
    @inbounds begin
        o = SV(h.origin[1],h.origin[2],h.origin[3])
        w = SV(h.widths[1],h.widths[2],h.widths[3])
        (o,
        o.+SV(   z,w[2],   z),
        o.+SV(w[1],w[2],   z),
        o.+SV(w[1],   z,   z),
        o.+SV(   z,   z,w[3]),
        o.+SV(   z,w[2],w[3]),
        o.+w,
        o.+SV(w[1],   z,w[3]))
    end
end

function interpolate_mt(c)
    @inbounds begin
        (sum(c)*0.125, # center
        (c[1]+c[2]+c[3]+c[4])*0.25,
        (c[1]+c[4]+c[5]+c[8])*0.25,
        (c[1]+c[2]+c[5]+c[6])*0.25,
        (c[5]+c[6]+c[7]+c[8])*0.25,
        (c[2]+c[3]+c[6]+c[7])*0.25,
        (c[3]+c[4]+c[7]+c[8])*0.25)
    end
end

function face_center_vertices(h::HyperRectangle{N,T}) where{N,T}
    SV = SVector{3,T}
    z = zero(T)
    @inbounds begin
        o = SV(h.origin[1],h.origin[2],h.origin[3])
        w = SV(h.widths[1],h.widths[2],h.widths[3])
        hw = w.*0.5
        c = center(h)
        (SV(c[1],c[2],c[3]),
        o.+SV(hw[1],hw[2],    z),
        o.+SV(hw[1],    z,hw[3]),
        o.+SV(z    ,hw[2],hw[3]),
        o.+SV(hw[1],hw[2], w[3]),
        o.+SV(hw[1], w[2],hw[3]),
        o.+SV( w[1],hw[2],hw[3]))
    end
end

center(rect::HyperRectangle) = rect.origin + 0.5 * rect.widths

function octsplit(h::HyperRectangle)
    ET = eltype(h.origin)
    VT = Vec{3,ET}
    z = zero(ET)
    hw = h.widths
    nw = h.widths.*ET(0.5)
    no = h.origin
    @inbounds begin
        (HyperRectangle(no, nw),
        HyperRectangle(no.+nw, nw),
        HyperRectangle(no.+VT(nw[1],z,z), nw),
        HyperRectangle(no.+VT(z,nw[2],z), nw),
        HyperRectangle(no.+VT(z,z,nw[3]), nw),
        HyperRectangle(no.+VT(nw[1],nw[2],z), nw),
        HyperRectangle(no.+VT(nw[1],z,nw[3]), nw),
        HyperRectangle(no.+VT(z,nw[2],nw[3]), nw))
    end
end

function get_iso_vals(f, val_store, points::NTuple{N,T}) where {N, T}
    @inbounds begin
        return ntuple(N) do i
            get_iso_vals(f,val_store,points[i])
        end
    end
end

function get_iso_vals(f, val_store, pt)
    # tuples are faster to hash than SVec
    key_tup = pt.data # tuple in SVectors/GeometryTypes, probably shoudl make more generic
    if haskey(val_store,key_tup)
        val_store[key_tup]
    else
        val_store[key_tup] = f(pt)::valtype(val_store)
    end
end

@inline function _get_interpindex(iso_vals, iso)
    cubeindex = iso_vals[1] < iso ? 0x01 : 0x00
    iso_vals[2] < iso && (cubeindex |= 0x02)
    iso_vals[3] < iso && (cubeindex |= 0x04)
    iso_vals[4] < iso && (cubeindex |= 0x08)
    iso_vals[5] < iso && (cubeindex |= 0x10)
    iso_vals[6] < iso && (cubeindex |= 0x20)
    iso_vals[7] < iso && (cubeindex |= 0x40)
    cubeindex
end

"""
    _get_cubeindex_pos(iso_vals, iso)

given `iso_vals` and iso, return an 8 bit value corresponding
to each corner of a cube. In each bit position,
0 indicates in the isosurface and 1 indicates outside the surface,
where the sign convention indicates positive inside the surface
"""
@inline function _get_interpindex_pos(iso_vals, iso)
    cubeindex = iso_vals[1] > iso ? 0x01 : 0x00
    iso_vals[2] > iso && (cubeindex |= 0x02)
    iso_vals[3] > iso && (cubeindex |= 0x04)
    iso_vals[4] > iso && (cubeindex |= 0x08)
    iso_vals[5] > iso && (cubeindex |= 0x10)
    iso_vals[6] > iso && (cubeindex |= 0x20)
    iso_vals[7] > iso && (cubeindex |= 0x40)
    cubeindex
end


function isosurface(f::Function, method::AdaptiveMarchingCubes, ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{3, Int};
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
        points = vertices(cell, VertType)

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
            vertlist = find_vertices_interp(points, iso_vals, cubeindex, method.iso, method.eps)

            # Create the triangle
            method.reduceverts && _mc_unique_triangles!(vts, fcs, vertlist, cubeindex, FaceType)
            !method.reduceverts && _mc_create_triangles!(vts, fcs, vertlist, cubeindex, FaceType)
        end
    end
    vts,fcs
end


@inline function vertPos(f, e, width, origin, ::Type{VertType}) where {VertType}
    T = eltype(eltype(VertType))

    ixs     = voxEdgeCrnrs[e]
    c1 = voxCrnrPos(VertType)[ixs[1]].*width .+ origin
    c2 = voxCrnrPos(VertType)[ixs[2]].*width .+ origin
    a = brent(f, c1, c2.-c1) # find root using brents method

    c1 + a.*(c2.-c1)
end


@inline function getVertId(f, e, width, vals, iso::Real, origin, vtsAry::Vector, vertex_store, eps::Real)

    VertType = eltype(vtsAry)

    # calculate vert position
    v = vertPos(f, e, width, origin, VertType)
    vt_key = v.data
    if haskey(vertex_store, vt_key)
        return vertex_store[vt_key]
    else
        push!(vtsAry, v)
        return vertex_store[vt_key] = length(vtsAry)
    end
end

slivers(x,y,z) = x == y || y == z || x == z

"""
    procVox(vals, iso::Real, x, y, z, nx, ny,
                    vts::Dict, vtsAry::Vector, fcs::Vector,
                    eps::Real)

Processes a voxel, adding any new vertices and faces to the given
containers as necessary.
"""
function procVox(f, vals, iso::Real, width, origin, vtsAry::Vector, vertex_store, fcs::Vector,
                 eps::Real, cubeindex)
    VertType = eltype(vtsAry)
    FaceType = eltype(fcs)
    # check each sub-tetrahedron in the voxel
    @inbounds for i = 1:6
        tIx = tetIx(i, cubeindex)
        (tIx == 0x00 || tIx == 0x0f) && continue

        e = tetTri[tIx]
        v1 = getVertId(f, voxEdgeId(i, e[1]), width, vals, iso, origin, vtsAry, vertex_store, eps)
        v2 = getVertId(f, voxEdgeId(i, e[2]), width, vals, iso, origin, vtsAry, vertex_store, eps)
        v3 = getVertId(f, voxEdgeId(i, e[3]), width, vals, iso, origin, vtsAry, vertex_store, eps)

        # add the face to the list
        !slivers(v1,v2,v3) && push!(fcs, FaceType(v1,v2,v3))

        # bail if there are no more faces
        iszero(e[4]) && continue
        v4 = getVertId(f, voxEdgeId(i, e[4]), width, vals, iso, origin, vtsAry, vertex_store, eps)
        v5 = getVertId(f, voxEdgeId(i, e[5]), width, vals, iso, origin, vtsAry, vertex_store, eps)
        v6 = getVertId(f, voxEdgeId(i, e[6]), width, vals, iso, origin, vtsAry, vertex_store, eps)
        !slivers(v4,v5,v6) && push!(fcs, FaceType(v4,v5,v6))

    end
end

function isosurface(f::Function, method::AdaptiveMarchingTetrahedra, ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{3, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2)) where {VertType, FaceType}

    ET = eltype(VertType)

    # arrays for vertices and faces
    vts = VertType[]
    fcs = FaceType[]

    # refinement queue
    # initialize with depth two
    refinement_queue = HyperRectangle{3,ET}[]
    for h in octsplit(HyperRectangle{3,ET}(origin,widths))
        append!(refinement_queue,octsplit(h))
    end

    val_store = Dict{NTuple{3,ET},ET}();
    vertex_store = Dict{NTuple{3,ET},ET}();

    @inbounds while true

        if isempty(refinement_queue)
            break
        end

        cell = pop!(refinement_queue)
        points = vertices_mt(cell, VertType)
        interp_points = face_center_vertices(cell)

        iso_vals = get_iso_vals(f,val_store,points)
        true_vals = get_iso_vals(f,val_store,interp_points)

        #Determine the index into the edge table which
        #tells us which vertices are inside of the surface
        cubeindex = method.insidepositive ? _get_cubeindex_pos(iso_vals, method.iso) : _get_cubeindex(iso_vals, method.iso)
        interpindex = method.insidepositive ? _get_interpindex_pos(true_vals, method.iso) : _get_interpindex(true_vals, method.iso)

        if (cubeindex == 0xff && interpindex == 0x7f) || (iszero(cubeindex) && iszero(interpindex))
            continue
        else
            value_interp = interpolate_mt(iso_vals)
            accurate = true
            for i = 1:7
                if !isapprox(true_vals[i], value_interp[i], rtol=method.rtol, atol=method.atol)
                    append!(refinement_queue, octsplit(cell))
                    accurate = false
                    break
                end
            end
            if accurate
                procVox(f, iso_vals, method.iso, cell.widths, cell.origin, vts, vertex_store, fcs, method.eps, cubeindex)
            end
        end
    end
    vts,fcs
end
