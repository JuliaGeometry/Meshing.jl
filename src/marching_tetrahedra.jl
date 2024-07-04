

#Look up Table
include("lut/mt.jl")

"""
    tetIx(tIx, cubeindex)

Determines which case in the triangle table we are dealing with
"""
function tetIx(tIx, cubeindex)
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
function vertId(e, x, y, z, nx, ny)
    dx, dy, dz = voxCrnrPos[voxEdgeCrnrs[e][1]]
    voxEdgeDir[e]+7*(x-0x01+dx+nx*(y-0x01+dy+ny*(z-0x01+dz)))
end

"""
    vertPos(e, x, y, z, vals::NTuple{8,T}, iso::Real, eps::Real, ::Type{VertType}) where {T<:Real, VertType}

Assuming an edge crossing, determines the point in space at which it
occurs.
eps represents the "bump" factor to keep vertices away from voxel
corners (thereby preventing degeneracies).
"""
function vertPos(e, x, y, z, xp, yp, zp, vals::V, iso, eps) where {V}
    T = eltype(vals)

    ixs     = voxEdgeCrnrs[e]
    srcVal  = vals[ixs[1]]
    tgtVal  = vals[ixs[2]]
    a       = min(max((iso-srcVal)/(tgtVal-srcVal), eps), one(T)-eps)
    b       = one(T)-a
    c1 = voxCrnrPos[ixs[1]]
    c2 = voxCrnrPos[ixs[2]]

    d = (xp[x+1], yp[y+1], zp[z+1]) .- (xp[x], yp[y], zp[z])
    (xp[x],yp[y],zp[z]) .+ (c1 .* b .+ c2 .* a) .* d
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
function getVertId(e, x, y, z, nx, ny,
                           vals, iso::Real,
                           xp, yp, zp,
                           vts::Dict,
                           vtsAry::Vector,
                           eps::Real)

    vId = vertId(e, x, y, z, nx, ny)
    haskey(vts, vId) && return vts[vId]

    # calculate vert position
    v = vertPos(e, x, y, z, xp, yp, zp, vals, iso, eps)
    push!(vtsAry, v)

    vts[vId] = length(vtsAry)

    return length(vtsAry)
end

"""
    voxEdgeId(subTetIx, tetEdgeIx)
Given a sub-tetrahedron case and a tetrahedron edge ID, determines the
corresponding voxel edge ID.
"""
function voxEdgeId(subTetIx, tetEdgeIx)
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
function procVox(vals, iso::Real, x, y, z, nx, ny, xp, yp, zp,
                 vts::Dict, vtsAry::Vector, fcs::Vector,
                 eps::Real, cubeindex)

    # check each sub-tetrahedron in the voxel
    @inbounds for i = 1:6
        tIx = tetIx(i, cubeindex)
        (tIx == 0x00 || tIx == 0x0f) && continue

        e = tetTri[tIx]

        # add the face to the list
        push!(fcs, (getVertId(voxEdgeId(i, e[1]), x, y, z, nx, ny, vals, iso, xp, yp, zp, vts, vtsAry, eps),
                    getVertId(voxEdgeId(i, e[2]), x, y, z, nx, ny, vals, iso, xp, yp, zp, vts, vtsAry, eps),
                    getVertId(voxEdgeId(i, e[3]), x, y, z, nx, ny, vals, iso, xp, yp, zp, vts, vtsAry, eps)))

        # bail if there are no more faces
        iszero(e[4]) && continue
        push!(fcs, (getVertId(voxEdgeId(i, e[4]), x, y, z, nx, ny, vals, iso, xp, yp, zp, vts, vtsAry, eps),
                    getVertId(voxEdgeId(i, e[5]), x, y, z, nx, ny, vals, iso, xp, yp, zp, vts, vtsAry, eps),
                    getVertId(voxEdgeId(i, e[6]), x, y, z, nx, ny, vals, iso, xp, yp, zp, vts, vtsAry, eps)))
    end
end

function isosurface(sdf::AbstractArray{T, 3}, method::MarchingTetrahedra, X=-1:1, Y=-1:1, Z=-1:1) where {T}

    FT = promote_type(eltype(first(X)), eltype(first(Y)), eltype(first(Z)), eltype(T), typeof(method.iso), typeof(method.eps))

    vts    = Dict{Int, Int}()
    fcs = NTuple{3,Int}[]
    vtsAry = NTuple{3,float(FT)}[]

    # process each voxel
    nx::Int, ny::Int, nz::Int = size(sdf)

    xp = LinRange(first(X), last(X), nx)
    yp = LinRange(first(Y), last(Y), ny)
    zp = LinRange(first(Z), last(Z), nz)

    @inbounds for i = 1:nx-1, j = 1:ny-1, k = 1:nz-1

        vals = (sdf[i,  j,  k  ],
                sdf[i,  j+1,k  ],
                sdf[i+1,j+1,k  ],
                sdf[i+1,j,  k  ],
                sdf[i,  j,  k+1],
                sdf[i,  j+1,k+1],
                sdf[i+1,j+1,k+1],
                sdf[i+1,j,  k+1])

        cubeindex = _get_cubeindex(vals, method.iso)

        _no_triangles(cubeindex) && continue

        procVox(vals, method.iso, i, j, k, nx, ny, xp, yp, zp, vts, vtsAry, fcs, method.eps, cubeindex)
    end

    vtsAry,fcs
end

function mt_vert_points(i, j, k, xp, yp, zp)
    ((xp[i], yp[j], zp[k]),
    (xp[i], yp[j+1], zp[k]),
    (xp[i+1], yp[j+1], zp[k]),
    (xp[i+1], yp[j], zp[k]),
    (xp[i], yp[j], zp[k+1]),
    (xp[i], yp[j+1], zp[k+1]),
    (xp[i+1], yp[j+1], zp[k+1]),
    (xp[i+1], yp[j], zp[k+1]))
end