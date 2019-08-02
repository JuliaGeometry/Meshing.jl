# API

## Quick Start

This package inherits the [mesh types](http://juliageometry.github.io/GeometryTypes.jl/latest/types.html#Meshes-1)
from [GeometryTypes.jl](https://github.com/JuliaGeometry/GeometryTypes.jl).

The algorithms operate on a `Function` or a `SignedDistanceField` and output a concrete `AbstractMesh`. For example:

```
using Meshing
using GeometryTypes
using LinearAlgebra: dot, norm
using FileIO

# Mesh an equation of sphere in the Axis-Aligned Bounding box starting
# at -1,-1,-1 and widths of 2,2,2
m = GLNormalMesh(HyperRectangle(Vec(-1,-1,-1.), Vec(2,2,2.)), MarchingCubes()) do v
    sqrt(sum(dot(v,v))) - 1
end

# save the Sphere as a PLY file
save("sphere.ply",m)
```

The general API is: ```(::Type{MT})(sdf::Function, method::AbstractMeshingAlgorithm) where {MT <: AbstractMesh}``` or ```(::Type{MT})(sdf::SignedDistanceField, method::AbstractMeshingAlgorithm) where {MT <: AbstractMesh}```


For a full listing of concrete `AbstractMesh` types see [GeometryTypes.jl mesh documentation](http://juliageometry.github.io/GeometryTypes.jl/latest/types.html#Meshes-1).

## Meshing Algorithms

Three meshing algorithms exist:
* `MarchingCubes()`
* `MarchingTetrahedra()`
* `NaiveSurfaceNets()`

Each takes an optional `iso` and `eps` parameter, e.g. `MarchingCubes(0.0,1e-6)`.

Here `iso` controls the offset for the boundary detection. By default this is set to 0. `eps` is the detection tolerance for a voxel edge intersection.

Below is a comparison of the algorithms:

| Algorithm          | Accurate | Manifold | Performance Penalty | Face Type |
|--------------------|----------|----------|---------------------|-----------|
| MarchingCubes      | Yes      | No       | ~2x                 | Triangle  |
| MarchingTetrahedra | Yes      | Yes      | ~3x                 | Triangle  |
| NaiveSurfaceNets   | No       | Yes      | 1x                  | Quad      |