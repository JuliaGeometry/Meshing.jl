
"""
    AbstractMeshingAlgorithm

Abstract type to specify an algorithm for isosurface extraction.
See:
- MarchingCubes
- MarchingTetrahedra
- NaiveSurfaceNets
"""
abstract type AbstractMeshingAlgorithm end

abstract type AbstractAdaptiveMeshingAlgorithm end


"""
    MarchingCubes(;iso=0.0)

Specifies the use of the Marching Cubes algorithm for isosurface extraction.
This algorithm provides a good balance between performance and vertex count.
In contrast to the other algorithms, vertices may be repeated, so mesh size
may be large and it will be difficult to extract topological/connectivity information.

- `iso` (default: 0.0) specifies the iso level to use for surface extraction.
"""
Base.@kwdef struct MarchingCubes{T} <: AbstractMeshingAlgorithm
    iso::T = 0.0
end

"""
    MarchingTetrahedra(;iso=0.0, eps=1e-3)

Specifies the use of the Marching Tetrahedra algorithm for isosurface extraction.
This algorithm has a roughly 2x performance penalty compared to Marching Cubes,
and generates more faces. However, each vertex is guaranteed to not be repeated,
making this algorithm useful for topological analysis.

- `iso` specifies the iso level to use for surface extraction.
- `eps` is the tolerence around a voxel corner to ensure manifold mesh generation.
"""
Base.@kwdef struct MarchingTetrahedra{T,E} <: AbstractMeshingAlgorithm
    iso::T = 0.0
    eps::E = 1e-3
end

"""
    MarchingCubes(iso=0.0, eps=1e-3)
    MarchingCubes(;iso=0.0, eps=1e-3)
    MarchingCubes(iso)
    MarchingCubes(iso,eps)

Specifies the use of the Marching Cubes algorithm for isosurface extraction.
This algorithm provides a good balance between performance and vertex count.
In contrast to the other algorithms, vertices may be repeated, so mesh size
may be large and it will be difficult to extract topological/connectivity information.

- `iso` (default: 0.0) specifies the iso level to use for surface extraction.
- `eps` (default: 1e-3) is the tolerence around a voxel corner to ensure manifold mesh generation.
- `reduceverts` (default: true) if true will merge vertices within a voxel to reduce mesh size by around 30% and with slight performance improvement.
"""
struct AdaptiveMarchingCubes{T,E} <: AbstractAdaptiveMeshingAlgorithm
    iso::T
    eps::T
    rtol::T
    atol::T
end


struct AdaptiveMarchingTetrahedra{T} <: AbstractAdaptiveMeshingAlgorithm
    iso::T
    eps::T
    rtol::T
    atol::T
end
