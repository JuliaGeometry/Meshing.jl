
struct MarchingCubes{T} <: AbstractMeshingAlgorithm
    iso::T
    eps::T
    reduceverts::Bool
end

MarchingCubes(;iso::T1=0.0, eps::T2=1e-3, reduceverts::Bool=true) where {T1, T2} = MarchingCubes{promote_type(T1, T2)}(iso, eps, reduceverts)
MarchingCubes(iso) = MarchingCubes(iso=iso)
MarchingCubes(iso,eps) = MarchingCubes(iso=iso,eps=eps)


struct MarchingTetrahedra{T} <: AbstractMeshingAlgorithm
    iso::T
    eps::T
    reduceverts::Bool
end

MarchingTetrahedra(;iso::T1=0.0, eps::T2=1e-3, reduceverts::Bool=true) where {T1, T2} = MarchingTetrahedra{promote_type(T1, T2)}(iso, eps, reduceverts)
MarchingTetrahedra(iso) = MarchingTetrahedra(iso=iso)
MarchingTetrahedra(iso,eps) = MarchingTetrahedra(iso=iso,eps=eps)


struct NaiveSurfaceNets{T} <: AbstractMeshingAlgorithm
    iso::T
    eps::T
    reduceverts::Bool
end

NaiveSurfaceNets(;iso::T1=0.0, eps::T2=1e-3, reduceverts::Bool=true) where {T1, T2} = NaiveSurfaceNets{promote_type(T1, T2)}(iso, eps, reduceverts)
NaiveSurfaceNets(iso) = NaiveSurfaceNets(iso=iso)
NaiveSurfaceNets(iso,eps) = NaiveSurfaceNets(iso=iso,eps=eps)
