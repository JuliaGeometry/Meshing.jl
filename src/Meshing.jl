module Meshing

using GeometryTypes,
      StaticArrays

"""
    _DEFAULT_SAMPLES = (24,24,24)

Global default sampling count for functions.
"""
const _DEFAULT_SAMPLES = (24,24,24)

include("algorithmtypes.jl")
include("geometrytypes_api.jl")
include("common.jl")
include("marching_tetrahedra.jl")
include("marching_cubes.jl")
include("surface_nets.jl")

export isosurface,
       MarchingCubes,
       MarchingTetrahedra,
       NaiveSurfaceNets

end # module
