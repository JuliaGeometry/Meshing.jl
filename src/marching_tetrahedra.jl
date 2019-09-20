

#Look up Table
include("lut/mt.jl")

"""
    tetIx(tIx, vals, iso::Real)

Determines which case in the triangle table we are dealing with
"""
@inline function tetIx(tIx, vals, iso::Real)
    v1 = vals[subTets[tIx][1]]
    v2 = vals[subTets[tIx][2]]
    v3 = vals[subTets[tIx][3]]
    v4 = vals[subTets[tIx][4]]
    ix = v1 < iso ? 0x01 : 0x00
    (v2 < iso) && (ix |= 0x02)
    (v3 < iso) && (ix |= 0x04)
    (v4 < iso) && (ix |= 0x08)
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
    voxEdgeDir[e]+7*(x-1+dx+nx*(y-1+dy+ny*(z-1+dz)))
end

"""
    vertPos(e, x, y, z, vals::NTuple{8,T}, iso::Real, eps::Real, ::Type{VertType}) where {T<:Real, VertType}

Assuming an edge crossing, determines the point in space at which it
occurs.
eps represents the "bump" factor to keep vertices away from voxel
corners (thereby preventing degeneracies).
"""
@inline function vertPos(e, x, y, z, vals::NTuple{8,T}, iso::Real, eps::Real, ::Type{VertType}) where {T<:Real, VertType}

    ixs     = voxEdgeCrnrs[e]
    srcVal  = vals[ixs[1]]
    tgtVal  = vals[ixs[2]]
    a       = min(max((iso-srcVal)/(tgtVal-srcVal), eps), one(T)-eps)
    b       = one(T)-a
    c1= voxCrnrPos(VertType)[ixs[1]]
    c2 = voxCrnrPos(VertType)[ixs[2]]

    VertType(x,y,z) + c1 .* b + c2.* a
end

"""
    getVertId(e, x, y, z, nx, ny,
                vals, iso::Real,
                vts::Dict,
                vtsAry::Vector,
                eps::Real, ::Type{VertType}) where {VertType}

Gets the vertex ID, adding it to the vertex dictionary if not already
present.
"""
@inline function getVertId(e, x, y, z, nx, ny,
                           vals, iso::Real,
                           vts::Dict,
                           vtsAry::Vector,
                           eps::Real, ::Type{VertType}) where {VertType}

    vId = vertId(e, x, y, z, nx, ny)
    # TODO we can probably immediately construct the vertex array here and use vert id to map to sequential ordering
    if !haskey(vts, vId)
        v = vertPos(e, x, y, z, vals, iso, eps, VertType)
        push!(vtsAry, v)
        vts[vId] = length(vtsAry)
    end
    vts[vId]
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
                    eps::Real,
                    ::Type{VertType}, ::Type{FaceType}) where {VertType, FaceType}

Processes a voxel, adding any new vertices and faces to the given
containers as necessary.
"""
function procVox(vals, iso::Real, x, y, z, nx, ny,
                 vts::Dict, vtsAry::Vector, fcs::Vector,
                 eps::Real,
                 ::Type{VertType}, ::Type{FaceType}) where {VertType, FaceType}

    # check each sub-tetrahedron in the voxel
    @inbounds for i = 1:6
        tIx = tetIx(i, vals, iso)
        (tIx == 0x00 || tIx == 0x0f) && continue

        e = tetTri[tIx]

        # add the face to the list
        push!(fcs, FaceType(
                    getVertId(voxEdgeId(i, e[1]), x, y, z, nx, ny, vals, iso, vts, vtsAry, eps, VertType),
                    getVertId(voxEdgeId(i, e[2]), x, y, z, nx, ny, vals, iso, vts, vtsAry, eps, VertType),
                    getVertId(voxEdgeId(i, e[3]), x, y, z, nx, ny, vals, iso, vts, vtsAry, eps, VertType)))

        # bail if there are no more faces
        e[4] == 0 && continue
        push!(fcs, FaceType(
                    getVertId(voxEdgeId(i, e[4]), x, y, z, nx, ny, vals, iso, vts, vtsAry, eps, VertType),
                    getVertId(voxEdgeId(i, e[5]), x, y, z, nx, ny, vals, iso, vts, vtsAry, eps, VertType),
                    getVertId(voxEdgeId(i, e[6]), x, y, z, nx, ny, vals, iso, vts, vtsAry, eps, VertType)))
    end
end


"""
Given a 3D array and an isovalue, extracts a mesh represention of the
an approximate isosurface by the method of marching tetrahedra.
"""
function isosurface(sdf::AbstractArray{T, 3}, method::MarchingTetrahedra, ::Type{VertType}=SVector{3,Float32}, ::Type{FaceType}=SVector{3, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2)) where {T, VertType, FaceType}
    vts        = Dict{Int, Int}()
    fcs        = FaceType[]
    vtsAry = VertType[]
    sizehint!(vts, div(length(sdf),8))
    sizehint!(vtsAry, div(length(sdf),8))
    sizehint!(fcs, div(length(sdf),4))
    # process each voxel
    nx::Int, ny::Int, nz::Int = size(sdf)
    @inbounds for k = 1:nz-1, j = 1:ny-1, i= 1:nx-1
        vals = (sdf[i, j, k],
                sdf[i, j+1, k],
                sdf[i+1, j+1, k],
                sdf[i+1, j, k],
                sdf[i, j, k+1],
                sdf[i, j+1, k+1],
                sdf[i+1, j+1, k+1],
                sdf[i+1, j, k+1])
        cubeindex = _get_cubeindex(vals,method.iso)
        if cubeindex != 0x00 && cubeindex != 0xff
            procVox(vals, method.iso, i, j, k, nx, ny, vts, vtsAry, fcs, method.eps, VertType, FaceType)
        end
    end

    _correct_vertices!(vtsAry, size(sdf), origin, widths, VertType)

    vtsAry,fcs
end

"""
    _correct_vertices!(vts, size, origin, widths, VertType)

The marchingTetrahedra function returns vertices on the (1-based) indices of the
SDF's data, ignoring its actual bounds. This function adjusts the vertices in
place so that they correspond to points within the SDF bounds.
"""
function _correct_vertices!(vts, size, origin, widths, ::Type{VertType}) where {VertType}
    s = widths ./ VertType(size .- 1)  # subtract 1 because an SDF with N points per side has N-1 cells
    for i in eachindex(vts)
        vts[i] = (vts[i] .- 1) .* s .+ origin  # subtract 1 to fix 1-indexing
    end
end

