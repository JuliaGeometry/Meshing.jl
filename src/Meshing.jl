module Meshing


"""
    _DEFAULT_SAMPLES = (24,24,24)

Global default sampling count for functions.
"""
const _DEFAULT_SAMPLES = (24,24,24)

include("algorithmtypes.jl")
include("common.jl")
include("marching_tetrahedra.jl")
include("marching_cubes.jl")
include("surface_nets.jl")
include("roots.jl")
#include("adaptive.jl")
#include("marching_tetrahedra_adaptive.jl")

export isosurface,
       MarchingCubes,
       MarchingTetrahedra,
       NaiveSurfaceNets

end # module
