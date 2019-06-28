module Meshing

using GeometryTypes
using AdaptiveDistanceFields
import RegionTrees
using StaticArrays

abstract type AbstractMeshingAlgorithm end

include("marching_tetrahedra.jl")
include("marching_cubes.jl")
include("surface_nets.jl")

export AdaptiveMarchingCubes,
       marching_cubes,
       MarchingCubes,
       MarchingTetrahedra,
       NaiveSurfaceNets

end # module
