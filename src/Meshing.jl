module Meshing

using StaticArrays

abstract type AbstractMeshingAlgorithm end

include("algorithmtypes.jl")
include("common.jl")
include("marching_tetrahedra.jl")
include("marching_cubes.jl")
include("surface_nets.jl")

export isosurface,
       MarchingCubes,
       MarchingTetrahedra,
       NaiveSurfaceNets

end # module
