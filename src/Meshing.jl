module Meshing

using GeometryTypes

abstract type AbstractMeshingAlgorithm end

include("marching_tetrahedra.jl")
include("marching_cubes.jl")
include("surface_nets.jl")

export marching_cubes,
       MarchingCubes,
       MarchingTetrahedra,
       NaiveSurfaceNets

end # module
