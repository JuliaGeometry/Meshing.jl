

#Look up Table
include("lut/mt.jl")

"""
Determines which case in the triangle table we are dealing with
"""
@inline function tetIx(tIx::IType, vals, iso::Real) where {IType <: Integer}
    v1 = vals[subTets[tIx][1]]
    v2 = vals[subTets[tIx][2]]
    v3 = vals[subTets[tIx][3]]
    v4 = vals[subTets[tIx][4]]
    ix = v1 < iso ? 1 : 0
    (v2 < iso) && (ix |= 2)
    (v3 < iso) && (ix |= 4)
    (v4 < iso) && (ix |= 8)
    ix
end

"""
Determines a unique integer ID associated with the edge. This is used
as a key in the vertex dictionary. It needs to be both unambiguous (no
two edges get the same index) and unique (every edge gets the same ID
regardless of which of its neighboring voxels is asking for it) in order
for vertex sharing to be implemented properly.
"""
@inline function vertId(e::IType, x::IType, y::IType, z::IType,
                nx::IType, ny::IType) where {IType <: Integer}
    dx, dy, dz = voxCrnrPos[voxEdgeCrnrs[e][1]]
    voxEdgeDir[e]+7*(x-1+dx+nx*(y-1+dy+ny*(z-1+dz)))
end

"""
Assuming an edge crossing, determines the point in space at which it
occurs.
eps represents the "bump" factor to keep vertices away from voxel
corners (thereby preventing degeneracies).
"""
@inline function vertPos(e::IType, x::IType, y::IType, z::IType,
                 vals::NTuple{8,T}, iso::Real, eps::Real) where {T<:Real, IType <: Integer}

    ixs     = voxEdgeCrnrs[e]
    srcVal  = vals[ixs[1]]
    tgtVal  = vals[ixs[2]]
    a       = min(max((iso-srcVal)/(tgtVal-srcVal), eps), one(T)-eps)
    b       = one(T)-a
    c1= voxCrnrPos[ixs[1]]
    c2 = voxCrnrPos[ixs[2]]

    Point{3,Float64}(x,y,z) + c1 .* b + c2.* a
end

"""
Gets the vertex ID, adding it to the vertex dictionary if not already
present.
"""
@inline function getVertId(e::IType, x::IType, y::IType, z::IType,
                   nx::IType, ny::IType,
                   vals, iso::Real,
                   vts::Dict{IType, Point{3,S}},
                   eps::Real) where {T <: Real, S <: Real, IType <: Integer}

    vId = vertId(e, x, y, z, nx, ny)
    # TODO we can probably immediately construct the vertex array here and use vert id to map to sequential ordering
    if !haskey(vts, vId)
        vts[vId] = vertPos(e, x, y, z, vals, iso, eps)
    end
    vId
end

"""
Given a sub-tetrahedron case and a tetrahedron edge ID, determines the
corresponding voxel edge ID.
"""
@inline function voxEdgeId(subTetIx, tetEdgeIx, IType)
    srcVoxCrnr = subTets[subTetIx][tetEdgeCrnrs[tetEdgeIx][1]]
    tgtVoxCrnr = subTets[subTetIx][tetEdgeCrnrs[tetEdgeIx][2]]
    return IType(voxEdgeIx[srcVoxCrnr][tgtVoxCrnr])
end

"""
Processes a voxel, adding any new vertices and faces to the given
containers as necessary.
"""
function procVox(vals, iso::Real,
                 x::IType, y::IType, z::IType,
                 nx::IType, ny::IType,
                 vts::Dict{IType, Point{3,S}}, fcs::Vector{Face{3,IType}},
                 eps::Real) where {T <: Real, S <: Real, IType <: Integer}

    # check each sub-tetrahedron in the voxel
    @inbounds for i::IType = 1:6
        tIx = tetIx(i, vals, iso)
        (tIx == 0 || tIx == 15) && continue

        e1 = tetTri[tIx][1]
        e2 = tetTri[tIx][2]
        e3 = tetTri[tIx][3]

        # add the face to the list
        push!(fcs, Face{3,IType}(
                    getVertId(voxEdgeId(i, e1, IType), x, y, z, nx, ny, vals, iso, vts, eps),
                    getVertId(voxEdgeId(i, e2, IType), x, y, z, nx, ny, vals, iso, vts, eps),
                    getVertId(voxEdgeId(i, e3, IType), x, y, z, nx, ny, vals, iso, vts, eps)))

        e1 = tetTri[tIx][4]
        # bail if there are no more faces
        e1 == 0 && continue
        e2 = tetTri[tIx][5]
        e3 = tetTri[tIx][6]
        push!(fcs, Face{3,IType}(
                    getVertId(voxEdgeId(i, e1, IType), x, y, z, nx, ny, vals, iso, vts, eps),
                    getVertId(voxEdgeId(i, e2, IType), x, y, z, nx, ny, vals, iso, vts, eps),
                    getVertId(voxEdgeId(i, e3, IType), x, y, z, nx, ny, vals, iso, vts, eps)))
    end
end


"""
Given a 3D array and an isovalue, extracts a mesh represention of the
an approximate isosurface by the method of marching tetrahedra.
"""
function marchingTetrahedra(lsf::AbstractArray{T,3}, iso::Real, eps::Real, indextype::Type{IT}) where {T<:Real, IT <: Integer}
    vertex_eltype = promote_type(T, typeof(iso), typeof(eps))
    vts        = Dict{indextype, Point{3,vertex_eltype}}()
    fcs        = Face{3,indextype}[]
    sizehint!(vts, div(length(lsf),8))
    sizehint!(fcs, div(length(lsf),4))
    # process each voxel
    (nx::indextype,ny::indextype,nz::indextype) = size(lsf)
    @inbounds for k::indextype = 1:nz-1, j::indextype = 1:ny-1, i::indextype = 1:nx-1
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
            procVox(vals, iso, i, j, k, nx, ny, vts, fcs, eps)
        end
    end

    (vts,fcs)
end

function isosurface(lsf, isoval, eps, indextype=Int, index_start=one(Int))
    # get marching tetrahedra version of the mesh
    (vts, fcs) = marchingTetrahedra(lsf, isoval, eps, indextype)
    # normalize the mesh representation
    vtD = Dict{indextype,indextype}()
    sizehint!(vtD, length(vts))
    k = index_start
    for x in keys(vts)
        vtD[x] = k
        k += one(indextype)
    end
    # rewrite the face array with the new vertex values
    for i in eachindex(fcs)
        f = fcs[i]
        fcs[i] = Face{3, indextype}(vtD[f[1]], vtD[f[2]], vtD[f[3]])
    end
    vtAry = collect(values(vts))

    (vtAry, fcs)
end

isosurface(lsf,isoval) = isosurface(lsf,isoval, convert(eltype(lsf), 0.001))

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

function (::Type{MT})(sdf::SignedDistanceField, method::MarchingTetrahedra) where {MT <: AbstractMesh}
    vts, fcs = isosurface(sdf.data, method.iso, method.eps)
    _correct_vertices!(vts, sdf)
    MT(vts, fcs)::MT
end

function (::Type{MT})(volume::Array{T, 3}, method::MarchingTetrahedra) where {MT <: AbstractMesh, T}
    vts, fcs = isosurface(volume, convert(T, method.iso), convert(T, method.eps))
    MT(vts, fcs)::MT
end
