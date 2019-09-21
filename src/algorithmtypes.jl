
"""
    AbstractMeshingAlgorithm

Abstract type to specify an algorithm for isosurface extraction.
See:
- MarchingCubes
- MarchingTetrahedra
- NaiveSurfaceNets
"""
abstract type AbstractMeshingAlgorithm end

"""
    MarchingCubes(iso=0.0, eps=1e-3, reduceverts=true, insidepositive=false)
    MarchingCubes(;iso=0.0, eps=1e-3, reduceverts=true, insidepositive=false)
    MarchingCubes(iso)
    MarchingCubes(iso,eps)

Specifies the use of the Marching Cubes algorithm for isosurface extraction.
This algorithm provides a good balance between performance and vertex count.
In contrast to the other algorithms, vertices may be repeated, so mesh size
may be large and it will be difficult to extract topological/connectivity information.

- `iso` (default: 0.0) specifies the iso level to use for surface extraction.
- `eps` (default: 1e-3) is the tolerence around a voxel corner to ensure manifold mesh generation.
- `reduceverts` (default: true) if true will merge vertices within a voxel to reduce mesh size by around 30% and with slight performance improvement.
- `insidepositive` (default: false) set true if the sign convention inside the surface is positive, common for NRRD and DICOM data
"""
struct MarchingCubes{T} <: AbstractMeshingAlgorithm
    iso::T
    eps::T
    reduceverts::Bool
    insidepositive::Bool
end

MarchingCubes(;iso::T1=0.0, eps::T2=1e-3, reduceverts::Bool=true, insidepositive::Bool=false) where {T1, T2} = MarchingCubes{promote_type(T1, T2)}(iso, eps, reduceverts, insidepositive)
MarchingCubes(iso) = MarchingCubes(iso=iso)
MarchingCubes(iso,eps) = MarchingCubes(iso=iso,eps=eps)


"""
    MarchingTetrahedra(iso=0.0, eps=1e-3, reduceverts=true, insidepositive=false)
    MarchingTetrahedra(;iso=0.0, eps=1e-3, reduceverts=true, insidepositive=false)
    MarchingTetrahedra(iso)
    MarchingTetrahedra(iso,eps)

Specifies the use of the Marching Tetrahedra algorithm for isosurface extraction.
This algorithm has a roughly 2x performance penalty compared to Marching Cubes,
and generates more faces. However, each vertex is guaranteed to not be repeated,
making this algorithm useful for topological analysis.

- `iso` specifies the iso level to use for surface extraction.
- `eps` is the tolerence around a voxel corner to ensure manifold mesh generation.
- `reduceverts` reserved for future use.
- `insidepositive` reserved for future use.
"""
struct MarchingTetrahedra{T} <: AbstractMeshingAlgorithm
    iso::T
    eps::T
    reduceverts::Bool
    insidepositive::Bool
end

MarchingTetrahedra(;iso::T1=0.0, eps::T2=1e-3, reduceverts::Bool=true, insidepositive::Bool=false) where {T1, T2} = MarchingTetrahedra{promote_type(T1, T2)}(iso, eps, reduceverts, insidepositive)
MarchingTetrahedra(iso) = MarchingTetrahedra(iso=iso)
MarchingTetrahedra(iso,eps) = MarchingTetrahedra(iso=iso,eps=eps)

"""
    NaiveSurfaceNets(iso=0.0, eps=1e-3, reduceverts=true, insidepositive=false)
    NaiveSurfaceNets(;iso=0.0, eps=1e-3, reduceverts=true, insidepositive=false)
    NaiveSurfaceNets(iso)
    NaiveSurfaceNets(iso,eps)

Specifies the use of the Naive Surface Nets algorithm for isosurface extraction.
This algorithm has a slight performance advantage in some cases over Marching Cubes.
Each vertex is guaranteed to not be repeated (useful for topological analysis),
however the algorithm does not guarantee accuracy and generates quad faces.

- `iso` specifies the iso level to use for surface extraction.
- `eps` is the tolerence around a voxel corner to ensure manifold mesh generation.
- `reduceverts` reserved for future use.
- `insidepositive` reserved for future use.
"""
struct NaiveSurfaceNets{T} <: AbstractMeshingAlgorithm
    iso::T
    eps::T
    reduceverts::Bool
    insidepositive::Bool
end

NaiveSurfaceNets(;iso::T1=0.0, eps::T2=1e-3, reduceverts::Bool=true, insidepositive::Bool=false) where {T1, T2} = NaiveSurfaceNets{promote_type(T1, T2)}(iso, eps, reduceverts, insidepositive)
NaiveSurfaceNets(iso) = NaiveSurfaceNets(iso=iso)
NaiveSurfaceNets(iso,eps) = NaiveSurfaceNets(iso=iso,eps=eps)


#
# Helper functions
#

default_face_length(::Union{MarchingCubes,MarchingTetrahedra}) = 3
default_face_length(::NaiveSurfaceNets) = 4
