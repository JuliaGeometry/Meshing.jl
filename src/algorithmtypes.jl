
"""
    AbstractMeshingAlgorithm

Abstract type to specify an algorithm for isosurface extraction.
See:
- [MarchingCubes](@ref)
- [MarchingTetrahedra](@ref)
"""
abstract type AbstractMeshingAlgorithm end


"""
    MarchingCubes(;iso=0.0, reduceverts=false)

Specifies the use of the Marching Cubes algorithm for isosurface extraction.
This algorithm provides a good balance between performance and vertex count.
In contrast to the other algorithms, vertices may be repeated, so mesh size
may be large and it will be difficult to extract topological/connectivity information.

- `iso` (default: 0.0) specifies the iso level to use for surface extraction.
- `reduceverts` (default: false) shares vertices between neighbouring voxels, giving a
  connected mesh with about a quarter of the vertices.
"""
Base.@kwdef struct MarchingCubes{T} <: AbstractMeshingAlgorithm
    iso::T = 0.0
    reduceverts::Bool = false
end
MarchingCubes(iso) = MarchingCubes(iso, false)  # the positional form from before `reduceverts`

"""
    MarchingTetrahedra(;iso=0.0, eps=1e-3)

Specifies the use of the Marching Tetrahedra algorithm for isosurface extraction.
This algorithm generates more faces. However, each vertex is guaranteed to not be repeated,
making this algorithm useful for topological analysis.

- `iso` specifies the iso level to use for surface extraction.
- `eps` is the tolerence around a voxel corner to ensure manifold mesh generation.
"""
Base.@kwdef struct MarchingTetrahedra{T,E} <: AbstractMeshingAlgorithm
    iso::T = 0.0
    eps::E = 1e-3
end

