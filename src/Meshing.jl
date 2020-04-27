module Meshing

using GeometryBasics,
      StaticArrays

"""
    _DEFAULT_SAMPLES = (24,24,24)

Global default sampling count for functions.
"""
const _DEFAULT_SAMPLES = (24,24,24)

include("distancefields.jl")
include("algorithmtypes.jl")
include("geometrybasics_api.jl")
include("common.jl")
include("marching_tetrahedra.jl")
include("marching_cubes.jl")
include("surface_nets.jl")

export isosurface,
       MarchingCubes,
       MarchingTetrahedra,
       NaiveSurfaceNets,
       SignedDistanceField

end # module
