
#
# General isosurface docstring
#

@doc """
    isosurface(V)
    isosurface(V, method::AbstractMeshingAlgorithm, X, Y, Z)

`isosurface` is the general interface to all isosurface extraction algorithms.

Returns: (Vector{NTuple{3, T}, }, Vector{NTuple{3, Int}})

Defaults:
`method` must be an instance of an `AbstractMeshingAlgorithm`, e.g.:
- MarchingCubes()
- MarchingTetrahedra()

If `isosurface` is called without a specified algorithm, it will default to MarchingCubes.

If a subtype of `AbstractArray` is specified, the mesh will by default be centered at the origin between
(-1,1) in each axis.


See also:
- [MarchingCubes](@ref)
- [MarchingTetrahedra](@ref)

"""
function isosurface(A, args...)
    isosurface(A, MarchingCubes(), args...)
end



