module Meshing

using GeometryTypes
using Roots
using ForwardDiff
using LinearAlgebra
using StaticArrays

abstract type AbstractMeshingAlgorithm end

include("marching_tetrahedra.jl")
include("marching_cubes.jl")
include("surface_nets.jl")
include("dual_contours.jl")

export marching_cubes,
       MarchingCubes,
       MarchingTetrahedra,
       NaiveSurfaceNets,
       DualContours

end # module
