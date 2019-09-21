
"""
    _get_cubeindex(iso_vals, iso)

given `iso_vals` and iso, return an 8 bit value corresponding
to each corner of a cube. In each bit position,
0 indicates in the isosurface and 1 indicates outside the surface,
where the sign convention indicates negative inside the surface
"""
@inline function _get_cubeindex(iso_vals, iso)
    cubeindex = iso_vals[1] < iso ? 0x01 : 0x00
    iso_vals[2] < iso && (cubeindex |= 0x02)
    iso_vals[3] < iso && (cubeindex |= 0x04)
    iso_vals[4] < iso && (cubeindex |= 0x08)
    iso_vals[5] < iso && (cubeindex |= 0x10)
    iso_vals[6] < iso && (cubeindex |= 0x20)
    iso_vals[7] < iso && (cubeindex |= 0x40)
    iso_vals[8] < iso && (cubeindex |= 0x80)
    cubeindex
end

"""
    _get_cubeindex_pos(iso_vals, iso)

given `iso_vals` and iso, return an 8 bit value corresponding
to each corner of a cube. In each bit position,
0 indicates in the isosurface and 1 indicates outside the surface,
where the sign convention indicates positive inside the surface
"""
@inline function _get_cubeindex_pos(iso_vals, iso)
    cubeindex = iso_vals[1] > iso ? 0x01 : 0x00
    iso_vals[2] > iso && (cubeindex |= 0x02)
    iso_vals[3] > iso && (cubeindex |= 0x04)
    iso_vals[4] > iso && (cubeindex |= 0x08)
    iso_vals[5] > iso && (cubeindex |= 0x10)
    iso_vals[6] > iso && (cubeindex |= 0x20)
    iso_vals[7] > iso && (cubeindex |= 0x40)
    iso_vals[8] > iso && (cubeindex |= 0x80)
    cubeindex
end

"""
    _determine_types(meshtype, fieldtype=Float64, facelen=3)

Given a subtype of AbstractMesh, determine the
type of vertex/point and face to use for internal computations.

Preference is given to the types specified by the Mesh call,
and will default to the `FieldType` for `SignedDistanceField`,
and Point{3,Float64}/Face{3,Int} for direct function sampling.
"""
function _determine_types(meshtype, fieldtype=Float64, facelen=3)
    # determine the point and face types
    # preference is given to the Mesh types
    # followed by SDF if unspecified
    if vertextype(meshtype) !== Any
        VertType = vertextype(meshtype)
    else
        VertType = Point{3, fieldtype}
    end
    if facetype(meshtype) !== Any
        FaceType = facetype(meshtype)
    else
        FaceType = Face{facelen, Int}
    end
    VertType, FaceType
end

#
# General isosurface docstring
#

@doc """

    function isosurface(sdf::AbstractArray{T, 3}, method::AbstractMeshingAlgorithm,
                         [ VertType = SVector{3,Float64} ], [ FaceType} = SVector{3, Int} ] ;
                         origin = SVector(-1.0,-1.0,-1.0), widths = SVector(2.0,2.0,2.0))

    function isosurface(f::Function, method::AbstractMeshingAlgorithm,
                         [ VertType = SVector{3,Float64} ], [FaceType = SVector{3, Int} ] ;
                         origin = SVector(-1.0,-1.0,-1.0), widths = SVector(2.0,2.0,2.0)
                         samples=(24,24,24))`

`isosurface` is the general interface to all isosurface extraction algorithms.

Returns: (Vector{VertType}, Vector{FaceType})

Defaults:
- VertType = SVector{3,Float64}
- FaceType = SVector{3, Int} ] ;
- origin = SVector(-1.0,-1.0,-1.0)
- widths = SVector(2.0,2.0,2.0)
- samples=(24,24,24) (function sampling only)

`method` must be an instance of an `AbstractMeshingAlgorithm`

If a subtype of `AbstractArray` is specified, the mesh will be default be centered at the origin between
(-1,1) in each axis. This may be overridden by specifying a new origin and widths for the axis-aligned bounding box
using keywords of the same names. For example if we want our vertices in the range of (0,1), we can specify `origin=SVector(0,0,0)`
and `widths = SVector(1,1,1)`.

If a function is specified, it will be uniformly sampled in each axis by the amount specified in `samples`.
The function is called the a single argument of `VertType`.

Performance Tips:
- ensure `VertType`, `origin`, and `widths` are all of the same type
- ensure the element type of `VertType` is the same as the specified isolevel

See also:
- MarchingCubes
- MarchingTetrahedra
- NaiveSurfaceNets

""" isosurface
