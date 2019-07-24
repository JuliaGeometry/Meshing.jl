

#Look up Table
include("lut/mt.jl")

"""
Determines which case in the triangle table we are dealing with
"""
function tetIx(tIx::IType, vals, iso::Real, vxidx::VoxelIndices{IType}) where {IType <: Integer}
    @inbounds v1 = vals[vxidx.subTets[tIx][1]]
    @inbounds v2 = vals[vxidx.subTets[tIx][2]]
    @inbounds v3 = vals[vxidx.subTets[tIx][3]]
    @inbounds v4 = vals[vxidx.subTets[tIx][4]]
    ifelse(v1 < iso, 1, 0) +
    ifelse(v2 < iso, 2, 0) +
    ifelse(v3 < iso, 4, 0) +
    ifelse(v4 < iso, 8, 0) + 1
end

"""
Determines a unique integer ID associated with the edge. This is used
as a key in the vertex dictionary. It needs to be both unambiguous (no
two edges get the same index) and unique (every edge gets the same ID
regardless of which of its neighboring voxels is asking for it) in order
for vertex sharing to be implemented properly.
"""
function vertId(e::IType, x::IType, y::IType, z::IType,
nx::IType, ny::IType, vxidx::VoxelIndices{IType}) where {IType <: Integer}
    @inbounds dx, dy, dz = vxidx.voxCrnrPos[vxidx.voxEdgeCrnrs[e][1]]
    vxidx.voxEdgeDir[e]+7*(x-1+dx+nx*(y-1+dy+ny*(z-1+dz)))
end

"""
Assuming an edge crossing, determines the point in space at which it
occurs.
eps represents the "bump" factor to keep vertices away from voxel
corners (thereby preventing degeneracies).
"""
function vertPos(e::IType, x::IType, y::IType, z::IType,
vals::NTuple{8,T}, iso::Real, eps::Real, vxidx::VoxelIndices{IType}) where {T<:Real, IType <: Integer}

    @inbounds ixs     = vxidx.voxEdgeCrnrs[e]
    @inbounds srcVal  = vals[ixs[1]]
    @inbounds tgtVal  = vals[ixs[2]]
    a       = min(max((iso-srcVal)/(tgtVal-srcVal), eps), one(T)-eps)
    b       = one(T)-a
    @inbounds c1x,c1y,c1z = vxidx.voxCrnrPos[ixs[1]]
    @inbounds c2x,c2y,c2z = vxidx.voxCrnrPos[ixs[2]]

    Point(
          x+b*c1x+a*c2x,
          y+b*c1y+a*c2y,
          z+b*c1z+a*c2z
    )
end

"""
Gets the vertex ID, adding it to the vertex dictionary if not already
present.
"""
function getVertId(e::IType, x::IType, y::IType, z::IType,
 nx::IType, ny::IType,
 vals, iso::Real,
 vts::Dict{IType, Point{3,S}},
 eps::Real, vxidx::VoxelIndices{IType}) where {T <: Real, S <: Real, IType <: Integer}

    vId = vertId(e, x, y, z, nx, ny, vxidx)
    if !haskey(vts, vId)
        vts[vId] = vertPos(e, x, y, z, vals, iso, eps, vxidx)
    end
    vId
end

"""
Given a sub-tetrahedron case and a tetrahedron edge ID, determines the
corresponding voxel edge ID.
"""
function voxEdgeId(subTetIx::IType, tetEdgeIx::IType, vxidx::VoxelIndices{IType}) where IType <: Integer
    @inbounds srcVoxCrnr::IType = vxidx.subTets[subTetIx][vxidx.tetEdgeCrnrs[tetEdgeIx][1]]
    @inbounds tgtVoxCrnr::IType = vxidx.subTets[subTetIx][vxidx.tetEdgeCrnrs[tetEdgeIx][2]]
    @inbounds v = vxidx.voxEdgeIx[srcVoxCrnr][tgtVoxCrnr]
    v
end

"""
Processes a voxel, adding any new vertices and faces to the given
containers as necessary.
"""
function procVox(vals, iso::Real,
                 x::IType, y::IType, z::IType,
                 nx::IType, ny::IType,
                 vts::Dict{IType, Point{3,S}}, fcs::Vector{Face{3,IType}},
                 eps::Real, vxidx::VoxelIndices{IType}) where {T <: Real, S <: Real, IType <: Integer}

    # check each sub-tetrahedron in the voxel
    for i::IType = 1:6
        tIx = tetIx(i, vals, iso, vxidx)
        for j::IType in 1:3:4
            @inbounds e1 = vxidx.tetTri[tIx][j]
            # bail if there are no more faces
            e1 == 0 && break
            @inbounds e2 = vxidx.tetTri[tIx][j+1]
            @inbounds e3 = vxidx.tetTri[tIx][j+2]

            # add the face to the list
            push!(fcs, Face{3,IType}(
                      getVertId(voxEdgeId(i, e1, vxidx), x, y, z, nx, ny, vals, iso, vts, eps, vxidx),
                      getVertId(voxEdgeId(i, e2, vxidx), x, y, z, nx, ny, vals, iso, vts, eps, vxidx),
                      getVertId(voxEdgeId(i, e3, vxidx), x, y, z, nx, ny, vals, iso, vts, eps, vxidx)))
        end
    end
end


"""
Given a 3D array and an isovalue, extracts a mesh represention of the
an approximate isosurface by the method of marching tetrahedra.
"""
function marchingTetrahedra(lsf::AbstractArray{T,3}, iso::Real, eps::Real, indextype::Type{IT}) where {T<:Real, IT <: Integer}
    vertex_eltype = promote_type(T, typeof(iso), typeof(eps))
    vts        = Dict{indextype, Point{3,vertex_eltype}}()
    fcs        = Array{Face{3,indextype}}(undef, 0)
    sizehint!(vts, div(length(lsf),8))
    sizehint!(fcs, div(length(lsf),4))
    vxidx = VoxelIndices{indextype}()
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
            procVox(vals, iso, i, j, k, nx, ny, vts, fcs, eps, vxidx)
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
    fcAry = Face{3,indextype}[Face{3,indextype}(vtD[f[1]], vtD[f[2]], vtD[f[3]]) for f in fcs]
    vtAry = collect(values(vts))

    (vtAry, fcAry)
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
