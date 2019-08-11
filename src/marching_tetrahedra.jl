

#Look up Table
include("lut/mt.jl")

"""
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
Determines a unique integer ID associated with the edge. This is used
as a key in the vertex dictionary. It needs to be both unambiguous (no
two edges get the same index) and unique (every edge gets the same ID
regardless of which of its neighboring voxels is asking for it) in order
for vertex sharing to be implemented properly.
"""
@inline function vertId(e, x, y, z, nx, ny)
    dx, dy, dz = voxCrnrPos[voxEdgeCrnrs[e][1]]
    voxEdgeDir[e]+7*(x-1+dx+nx*(y-1+dy+ny*(z-1+dz)))
end

"""
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
    c1= voxCrnrPos[ixs[1]]
    c2 = voxCrnrPos[ixs[2]]

    VertType(x,y,z) + c1 .* b + c2.* a
end

"""
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
Given a sub-tetrahedron case and a tetrahedron edge ID, determines the
corresponding voxel edge ID.
"""
@inline function voxEdgeId(subTetIx, tetEdgeIx)
    srcVoxCrnr = subTets[subTetIx][tetEdgeCrnrs[tetEdgeIx][1]]
    tgtVoxCrnr = subTets[subTetIx][tetEdgeCrnrs[tetEdgeIx][2]]
    return voxEdgeIx[srcVoxCrnr][tgtVoxCrnr]
end

"""
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
function marchingTetrahedra(lsf::AbstractArray{T,3}, iso::Real, eps::Real, ::Type{VertType}, ::Type{FaceType}) where {T<:Real, VertType, FaceType}
    vts        = Dict{Int, Int}()
    fcs        = FaceType[]
    vtsAry = VertType[]
    sizehint!(vts, div(length(lsf),8))
    sizehint!(vtsAry, div(length(lsf),8))
    sizehint!(fcs, div(length(lsf),4))
    # process each voxel
    nx::Int, ny::Int, nz::Int = size(lsf)
    @inbounds for k = 1:nz-1, j = 1:ny-1, i= 1:nx-1
        vals = (lsf[i, j, k],
                lsf[i, j+1, k],
                lsf[i+1, j+1, k],
                lsf[i+1, j, k],
                lsf[i, j, k+1],
                lsf[i, j+1, k+1],
                lsf[i+1, j+1, k+1],
                lsf[i+1, j, k+1])
        cubeindex = _get_cubeindex(vals,iso)
        if cubeindex != 0x00 && cubeindex != 0xff
            procVox(vals, iso, i, j, k, nx, ny, vts, vtsAry, fcs, eps, VertType, FaceType)
        end
    end

    (vtsAry,fcs)
end

"""
The marchingTetrahedra function returns vertices on the (1-based) indices of the
SDF's data, ignoring its actual bounds. This function adjusts the vertices in
place so that they correspond to points within the SDF bounds.
"""
function _correct_vertices!(vts, sdf::SignedDistanceField)
    bounds = HyperRectangle(sdf)
    orig = origin(bounds)
    w = widths(bounds)
    s = w ./ Point(size(sdf) .- 1)  # subtract 1 because an SDF with N points per side has N-1 cells
    for i in eachindex(vts)
        vts[i] = (vts[i] .- 1) .* s .+ orig  # subtract 1 to fix 1-indexing
    end
end

struct MarchingTetrahedra{T} <: AbstractMeshingAlgorithm
    iso::T
    eps::T
    reduceverts::Bool
end

MarchingTetrahedra(;iso::T1=0.0, eps::T2=1e-3, reduceverts::Bool=true) where {T1, T2} = MarchingTetrahedra{promote_type(T1, T2)}(iso, eps, reduceverts)
MarchingTetrahedra(iso) = MarchingTetrahedra(iso=iso)
MarchingTetrahedra(iso,eps) = MarchingTetrahedra(iso=iso,eps=eps)

function (::Type{MT})(sdf::SignedDistanceField{3,ST,FT}, method::MarchingTetrahedra) where {ST, FT, MT <: AbstractMesh}
    vertex_eltype = promote_type(FT, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype)
    vts, fcs = marchingTetrahedra(sdf.data, method.iso, method.eps, VertType, FaceType)
    _correct_vertices!(vts, sdf)
    MT(vts, fcs)::MT
end

function (::Type{MT})(volume::Array{T, 3}, method::MarchingTetrahedra) where {MT <: AbstractMesh, T}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype)
    vts, fcs = marchingTetrahedra(volume, convert(T, method.iso), convert(T, method.eps), VertType, FaceType)
    MT(vts, fcs)::MT
end
