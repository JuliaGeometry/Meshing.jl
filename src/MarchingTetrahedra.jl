"""
Marching Tetrahedra is an algorithm for extracting a triangular
mesh representation of an isosurface of a scalar volumetric
function sampled on a rectangular grid.

We divide the cube into six tetrahedra. [It is possible to divide
a cube into five tetrahedra, but not in a way that a translated
version of the division would share face diagonals. (It reqires a
reflection.)]
"""
module MarchingTetrahedra

using GeometryTypes

"""
Voxel corner and edge indexing conventions

        Z
        |

        5------5------6              Extra edges not drawn
       /|            /|              -----------
      8 |           6 |              - face diagonals
     /  9          /  10                - 13: 1 to 3
    8------7------7   |                 - 14: 1 to 8
    |   |         |   |                 - 15: 1 to 6
    |   1------1--|---2  -- Y           - 16: 5 to 7
    12 /          11 /                  - 17: 2 to 7
    | 4           | 2                   - 18: 4 to 7
    |/            |/                 - body diagonal
    4------3------3                     - 19: 1 to 7

  /
 X
"""

immutable VoxelIndices{T <: Integer}
    voxCrnrPos::NTuple{8,Vec{3,T}}
    voxEdgeCrnrs::NTuple{19, Vec{2,T}}
    voxEdgeDir::NTuple{19,T}
    voxEdgeIx::Matrix{T}
    subTets::Matrix{T}
    tetEdgeCrnrs::Matrix{T}
    tetTri::Matrix{T}

    function VoxelIndices()
        VT3 = Vec{3,T}
        VT2 = Vec{2,T}
        voxCrnrPos = (VT3(0, 0, 0),
                      VT3(0, 1, 0),
                      VT3(1, 1, 0),
                      VT3(1, 0, 0),
                      VT3(0, 0, 1),
                      VT3(0, 1, 1),
                      VT3(1, 1, 1),
                      VT3(1, 0, 1))
        # the voxel IDs at either end of the tetrahedra edges, by edge ID
        voxEdgeCrnrs = (VT2(1, 2),
                        VT2(2, 3),
                        VT2(4, 3),
                        VT2(1, 4),
                        VT2(5, 6),
                        VT2(6, 7),
                        VT2(8, 7),
                        VT2(5, 8),
                        VT2(1, 5),
                        VT2(2, 6),
                        VT2(3, 7),
                        VT2(4, 8),
                        VT2(1, 3),
                        VT2(1, 8),
                        VT2(1, 6),
                        VT2(5, 7),
                        VT2(2, 7),
                        VT2(4, 7),
                        VT2(1, 7))

        # direction codes:
        # 0 => +x, 1 => +y, 2 => +z,
        # 3 => +xy, 4 => +xz, 5 => +yz, 6 => +xyz
        voxEdgeDir = convert(NTuple{19,T}, (1,0,1,0,1,0,1,0,2,2,2,2,3,4,5,3,4,5,6))

        # For a pair of corner IDs, the edge ID joining them
        # 0 denotes a pair with no edge
        voxEdgeIx = [0  1  13 4  9  15 19 14;
                     1  0  2  0  0  10 17 0;
                     13 2  0  3  0  0  11 0;
                     4  0  3  0  0  0  18 12;
                     9  0  0  0  0  5  16 8;
                     15 10 0  0  5  0  6  0;
                     19 17 11 18 16 6  0  7;
                     14 0  0  12 8  0  7  0]

        # voxel corners that comprise each of the six tetrahedra
        subTets = [   1 3 2 7;
                      1 8 4 7;
                      1 4 3 7;
                      1 2 6 7;
                      1 5 8 7;
                      1 6 5 7]'
        # tetrahedron corners for each edge (indices 1-4)
        tetEdgeCrnrs = [1 2;
                        2 3;
                        1 3;
                        1 4;
                        2 4;
                        3 4]'

        # triangle cases for a given tetrahedron edge code
        tetTri = T[ 0 0 0 0 0 0;
                    1 3 4 0 0 0;
                    1 5 2 0 0 0;
                    3 5 2 3 4 5;
                    2 6 3 0 0 0;
                    1 6 4 1 2 6;
                    1 5 6 1 6 3;
                    4 5 6 0 0 0;
                    4 6 5 0 0 0;
                    1 6 5 1 3 6;
                    1 4 6 1 6 2;
                    2 3 6 0 0 0;
                    3 2 5 3 5 4;
                    1 2 5 0 0 0;
                    1 4 3 0 0 0;
                    0 0 0 0 0 0]'

        new(voxCrnrPos,
            voxEdgeCrnrs,
            voxEdgeDir,
            voxEdgeIx,
            subTets,
            tetEdgeCrnrs,
            tetTri)
    end
end
# (X,Y,Z)-coordinates for each voxel corner ID


# Checks if a voxel has faces. Should be false for most voxels.
# This function should be made as fast as possible.
function hasFaces{T<:Real}(vals::Vector{T}, iso::T)
    if vals[1] < iso
        for i=2:8
            vals[i] >= iso && return true
        end
    else
        for i=2:8
            vals[i] <  iso && return true
        end
    end
    false
end

# Determines which case in the triangle table we are dealing with
function tetIx{T<:Real, IType <: Integer}(tIx::IType, vals::Vector{T}, iso::T, vxidx::VoxelIndices{IType})
    ifelse(vals[vxidx.subTets[1, tIx]] < iso, 1, 0) +
    ifelse(vals[vxidx.subTets[2, tIx]] < iso, 2, 0) +
    ifelse(vals[vxidx.subTets[3, tIx]] < iso, 4, 0) +
    ifelse(vals[vxidx.subTets[4, tIx]] < iso, 8, 0) + 1
end

# Determines a unique integer ID associated with the edge. This is used
# as a key in the vertex dictionary. It needs to be both unambiguous (no
# two edges get the same index) and unique (every edge gets the same ID
# regardless of which of its neighboring voxels is asking for it) in order
# for vertex sharing to be implemented properly.
function vertId{IType <: Integer}(e::IType, x::IType, y::IType, z::IType,
                nx::IType, ny::IType, vxidx::VoxelIndices{IType})
    dx = vxidx.voxCrnrPos[vxidx.voxEdgeCrnrs[e][1]]
    vxidx.voxEdgeDir[e]+7*(x-1+dx[1]+nx*(y-1+dx[2]+ny*(z-1+dx[3])))
end

# Assuming an edge crossing, determines the point in space at which it
# occurs.
# eps represents the "bump" factor to keep vertices away from voxel
# corners (thereby preventing degeneracies).
function vertPos{T<:Real, IType <: Integer}(e::IType, x::IType, y::IType, z::IType,
                          vals::Vector{T}, iso::T, eps::T, vxidx::VoxelIndices{IType})

    ixs     = vxidx.voxEdgeCrnrs[e]
    srcVal  = vals[ixs[1]]
    tgtVal  = vals[ixs[2]]
    a       = (iso-srcVal)/(tgtVal-srcVal)
    a       = min(max(a, eps), one(T)-eps)
    b       = one(T)-a
    corner1 = vxidx.voxCrnrPos[ixs[1]]
    corner2 = vxidx.voxCrnrPos[ixs[2]]

    Point{3,T}(
          x+b*corner1[1]+a*corner2[1],
          y+b*corner1[2]+a*corner2[2],
          z+b*corner1[3]+a*corner2[3]
    )
end

# Gets the vertex ID, adding it to the vertex dictionary if not already
# present.
function getVertId{T<:Real, IType <: Integer}(e::IType, x::IType, y::IType, z::IType,
                            nx::IType, ny::IType,
                            vals::Vector{T}, iso::T,
                            vts::Dict{IType, Point{3,T}},
                            eps::T, vxidx::VoxelIndices{IType})

    vId = vertId(e, x, y, z, nx, ny, vxidx)
    if !haskey(vts, vId)
        vts[vId] = vertPos(e, x, y, z, vals, iso, eps, vxidx)
    end
    vId
end

# Given a sub-tetrahedron case and a tetrahedron edge ID, determines the
# corresponding voxel edge ID.
function voxEdgeId{IType <: Integer}(subTetIx::IType, tetEdgeIx::IType, vxidx::VoxelIndices{IType})
    srcVoxCrnr::IType = vxidx.subTets[vxidx.tetEdgeCrnrs[1, tetEdgeIx], subTetIx]
    tgtVoxCrnr::IType = vxidx.subTets[vxidx.tetEdgeCrnrs[2, tetEdgeIx], subTetIx]
    vxidx.voxEdgeIx[srcVoxCrnr, tgtVoxCrnr]
end

# Processes a voxel, adding any new vertices and faces to the given
# containers as necessary.
function procVox{T<:Real, IType <: Integer}(vals::Vector{T}, iso::T,
                          x::IType, y::IType, z::IType,
                          nx::IType, ny::IType,
                          vts::Dict{IType, Point{3,T}}, fcs::Vector{Face{3,IType,0}},
                          eps::T, vxidx::VoxelIndices{IType})

    # check each sub-tetrahedron in the voxel
    for i::IType = 1:6
        tIx = tetIx(i, vals, iso, vxidx)
        for j::IType in 1:3:4
            e1 = vxidx.tetTri[j, tIx]
            # bail if there are no more faces
            e1 == 0 && break
            e2 = vxidx.tetTri[j+1,tIx]
            e3 = vxidx.tetTri[j+2,tIx]

            # add the face to the list
            push!(fcs, Face{3,IType,0}(
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
function marchingTetrahedra{T<:Real, IT <: Integer}(lsf::AbstractArray{T,3}, iso::T, eps::T, indextype::Type{IT})
    vts        = Dict{indextype, Point{3,T}}()
    fcs        = Array(Face{3,indextype,0}, 0)
    sizehint!(vts, div(length(lsf),8))
    sizehint!(fcs, div(length(lsf),4))
    const vxidx = VoxelIndices{indextype}()
    # process each voxel
    (nx::indextype,ny::indextype,nz::indextype) = size(lsf)
    vals = zeros(T, 8)
    @inbounds for k::indextype = 1:nz-1, j::indextype = 1:ny-1, i::indextype = 1:nx-1
        for l::indextype=1:8
            vals[l] = lsf[i+vxidx.voxCrnrPos[l][1], j+vxidx.voxCrnrPos[l][2], k+vxidx.voxCrnrPos[l][3]]
        end
        if hasFaces(vals,iso)
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
    fcAry = Face{3,indextype, index_start-1}[Face{3,indextype, index_start-1}(vtD[f[1]], vtD[f[2]], vtD[f[3]]) for f in fcs]
    vtAry = collect(values(vts))

    (vtAry, fcAry)
end

isosurface(lsf,isoval) = isosurface(lsf,isoval, convert(eltype(lsf), 0.001))


function call{MT <: AbstractMesh, T}(::Type{MT}, volume::Array{T, 3}, iso_val::Real, eps_val=0.001)
    iso_val = convert(T, iso_val)
    eps_val = convert(T, eps_val)
    vts, fcs = isosurface(volume, iso_val, eps_val)
    MT(vts, fcs)
end

end # module
