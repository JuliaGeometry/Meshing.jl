

#Look up Table
include("lut/mt.jl")

"""
    tetIx(tIx, cubeindex)

Determines which case in the triangle table we are dealing with
"""
@inline function tetIx(tIx, cubeindex)
    @inbounds v1 = subTetsMask[tIx][1]
    @inbounds v2 = subTetsMask[tIx][2]
    ix = 0x01 & cubeindex | (0x40 & cubeindex) >> 0x03
    !iszero(v1 & cubeindex) && (ix |= 0x02)
    !iszero(v2 & cubeindex) && (ix |= 0x04)
    ix
end

"""
    vertId(e, x, y, z, nx, ny)

Determines a unique integer ID associated with the edge. This is used
as a key in the vertex dictionary. It needs to be both unambiguous (no
two edges get the same index) and unique (every edge gets the same ID
regardless of which of its neighboring voxels is asking for it) in order
for vertex sharing to be implemented properly.
"""
@inline function vertId(e, x, y, z, nx, ny)
    dx, dy, dz = voxCrnrPosInt[voxEdgeCrnrs[e][1]]
    voxEdgeDir[e]+7*(x-0x01+dx+nx*(y-0x01+dy+ny*(z-0x01+dz)))
end

"""
    vertPos(e, x, y, z, vals::NTuple{8,T}, iso::Real, eps::Real, ::Type{VertType}) where {T<:Real, VertType}

Assuming an edge crossing, determines the point in space at which it
occurs.
eps represents the "bump" factor to keep vertices away from voxel
corners (thereby preventing degeneracies).
"""
@inline function vertPos(e, x, y, z, scale, origin, vals::V, iso, eps, ::Type{VertType}) where {V, VertType}
    T = eltype(vals)

    ixs     = voxEdgeCrnrs[e]
    srcVal  = vals[ixs[1]]
    tgtVal  = vals[ixs[2]]
    a       = min(max((iso-srcVal)/(tgtVal-srcVal), eps), one(T)-eps)
    b       = one(T)-a
    c1 = voxCrnrPos(VertType)[ixs[1]]
    c2 = voxCrnrPos(VertType)[ixs[2]]

    ((VertType(x,y,z) + c1 .* b + c2.* a) .- 1) .* scale .+ origin
end

"""
    getVertId(e, x, y, z, nx, ny,
                vals, iso::Real,
                vts::Dict,
                vtsAry::Vector,
                eps::Real)

Gets the vertex ID, adding it to the vertex dictionary if not already
present.
"""
@inline function getVertId(e, x, y, z, nx, ny,
                           vals, iso::Real,
                           scale, origin,
                           vts::Dict,
                           vtsAry::Vector,
                           eps::Real,
                           reduceverts::Bool)

    VertType = eltype(vtsAry)
    if reduceverts
        vId = vertId(e, x, y, z, nx, ny)
        haskey(vts, vId) && return vts[vId]
    end

    # calculate vert position
    v = vertPos(e, x, y, z, scale, origin, vals, iso, eps, VertType)
    push!(vtsAry, v)

    # if deduplicting, push to dict
    if reduceverts
        vts[vId] = length(vtsAry)
    end

    return length(vtsAry)
end

"""
    voxEdgeId(subTetIx, tetEdgeIx)
Given a sub-tetrahedron case and a tetrahedron edge ID, determines the
corresponding voxel edge ID.
"""
@inline function voxEdgeId(subTetIx, tetEdgeIx)
    srcVoxCrnr = subTets[subTetIx][tetEdgeCrnrs[tetEdgeIx][1]]
    tgtVoxCrnr = subTets[subTetIx][tetEdgeCrnrs[tetEdgeIx][2]]
    return voxEdgeIx[srcVoxCrnr][tgtVoxCrnr]
end

"""
    procVox(vals, iso::Real, x, y, z, nx, ny,
                    vts::Dict, vtsAry::Vector, fcs::Vector,
                    eps::Real)

Processes a voxel, adding any new vertices and faces to the given
containers as necessary.
"""
function procVox(vals, iso::Real, x, y, z, nx, ny, scale, origin,
                 vts::Dict, vtsAry::Vector, fcs::Vector,
                 eps::Real, cubeindex, reduceverts)
    VertType = eltype(vtsAry)
    FaceType = eltype(fcs)
    # check each sub-tetrahedron in the voxel
    @inbounds for i = 1:6
        tIx = tetIx(i, cubeindex)
        (tIx == 0x00 || tIx == 0x0f) && continue

        e = tetTri[tIx]

        # add the face to the list
        push!(fcs, FaceType(
                    getVertId(voxEdgeId(i, e[1]), x, y, z, nx, ny, vals, iso, scale, origin, vts, vtsAry, eps, reduceverts),
                    getVertId(voxEdgeId(i, e[2]), x, y, z, nx, ny, vals, iso, scale, origin, vts, vtsAry, eps, reduceverts),
                    getVertId(voxEdgeId(i, e[3]), x, y, z, nx, ny, vals, iso, scale, origin, vts, vtsAry, eps, reduceverts)))

        # bail if there are no more faces
        iszero(e[4]) && continue
        push!(fcs, FaceType(
                    getVertId(voxEdgeId(i, e[4]), x, y, z, nx, ny, vals, iso, scale, origin, vts, vtsAry, eps, reduceverts),
                    getVertId(voxEdgeId(i, e[5]), x, y, z, nx, ny, vals, iso, scale, origin, vts, vtsAry, eps, reduceverts),
                    getVertId(voxEdgeId(i, e[6]), x, y, z, nx, ny, vals, iso, scale, origin, vts, vtsAry, eps, reduceverts)))
    end
end

function isosurface(sdf::AbstractArray{T, 3}, method::MarchingTetrahedra, ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{3, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2)) where {T, VertType, FaceType}

    vts    = Dict{Int, Int}()
    fcs    = FaceType[]
    vtsAry = VertType[]
    sizehint!(vts, div(length(sdf),8))
    sizehint!(vtsAry, div(length(sdf),8))
    sizehint!(fcs, div(length(sdf),4))
    # process each voxel
    scale = widths ./ VertType(size(sdf) .- 1)
    nx::Int, ny::Int, nz::Int = size(sdf)

    @inbounds for k = 1:nz-1, j = 1:ny-1, i = 1:nx-1

        vals = (sdf[i,  j,  k  ],
                sdf[i,  j+1,k  ],
                sdf[i+1,j+1,k  ],
                sdf[i+1,j,  k  ],
                sdf[i,  j,  k+1],
                sdf[i,  j+1,k+1],
                sdf[i+1,j+1,k+1],
                sdf[i+1,j,  k+1])

        cubeindex = method.insidepositive ? _get_cubeindex_pos(vals, method.iso) : _get_cubeindex(vals, method.iso)

        _no_triangles(cubeindex) && continue

        procVox(vals, method.iso, i, j, k, nx, ny, scale, origin, vts, vtsAry, fcs, method.eps, cubeindex, method.reduceverts)
    end

    vtsAry,fcs
end

function isosurface(f::Function, method::MarchingTetrahedra,
                    ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{3, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2),
                    samples::NTuple{3,T}=_DEFAULT_SAMPLES) where {T, VertType, FaceType}

    vts    = Dict{Int, Int}()
    fcs    = FaceType[]
    vtsAry = VertType[]
    sizehint!(vts, div(prod(samples),8))
    sizehint!(vtsAry, div(prod(samples),8))
    sizehint!(fcs, div(prod(samples),4))
    # process each voxel
    scale = widths ./ (VertType(samples...) .-1)
    nx::Int, ny::Int, nz::Int = samples
    zv = zero(eltype(VertType))
    vals = (zv,zv,zv,zv,zv,zv,zv,zv)

    @inbounds for k = 1:nz, j = 1:ny, i = 1:nx
        points = (VertType(i-1,j-1,k-1).* scale .+ origin,
                  VertType(i-1,j  ,k-1).* scale .+ origin,
                  VertType(i  ,j  ,k-1).* scale .+ origin,
                  VertType(i  ,j-1,k-1).* scale .+ origin,
                  VertType(i-1,j-1,k  ).* scale .+ origin,
                  VertType(i-1,j  ,k  ).* scale .+ origin,
                  VertType(i  ,j  ,k  ).* scale .+ origin,
                  VertType(i  ,j-1,k  ).* scale .+ origin)
        if i == 1
            vals = (f(points[1]),
                    f(points[2]),
                    f(points[3]),
                    f(points[4]),
                    f(points[5]),
                    f(points[6]),
                    f(points[7]),
                    f(points[8]))
        else
            vals = (vals[4],
                    vals[3],
                    f(points[3]),
                    f(points[4]),
                    vals[8],
                    vals[7],
                    f(points[7]),
                    f(points[8]))
        end

        cubeindex = method.insidepositive ? _get_cubeindex_pos(vals, method.iso) : _get_cubeindex(vals, method.iso)

        _no_triangles(cubeindex) && continue

        procVox(vals, method.iso, i, j, k, nx, ny, scale, origin, vts, vtsAry, fcs, method.eps, cubeindex, method.reduceverts)
    end

    vtsAry,fcs
end
