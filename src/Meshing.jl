VERSION >= v"0.4.0-dev+6521" && __precompile__()

module Meshing

using GeometryTypes

abstract type AbstractMeshingAlgorithm end

include("marching_tetrahedra.jl")
include("marching_cubes.jl")

export marching_cubes,
       MarchingCubes,
       MarchingTetrahedra

end # module
