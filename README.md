# Meshing

[![Build Status](https://travis-ci.org/JuliaGeometry/Meshing.jl.svg)](https://travis-ci.org/JuliaGeometry/Meshing.jl)
[![codecov.io](http://codecov.io/github/JuliaGeometry/Meshing.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaGeometry/Meshing.jl?branch=master)

This package provides meshing algorithms for use on distance fields.

Including:
* [Marching Tetrahedra](https://en.wikipedia.org/wiki/Marching_tetrahedra)
* [Marching Cubes](https://en.wikipedia.org/wiki/Marching_cubes)
* [Naive Surface Nets](https://0fps.net/2012/07/12/smooth-voxel-terrain-part-2/)

## Interface

This package is tightly integrated with [GeometryTypes.jl](https://github.com/JuliaGeometry/GeometryTypes.jl).

All algorithms operate on `SignedDistanceField` and output a concrete `AbstractMesh`. For example:

```
using Meshing
using GeometryTypes
using LinearAlgebra: dot, norm
using FileIO

# generate an SDF of a sphere
sdf_sphere = SignedDistanceField(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.))) do v
    sqrt(sum(dot(v,v))) - 1 # sphere
end

m = GLNormalMesh(sdf_sphere, MarchingCubes())

save("sphere.ply",m)
```

The general API is ``(::Type{MT})(sdf::SignedDistanceField, method::AbstractMeshingAlgorithm) where {MT <: AbstractMesh}``

For a full listing of concrete `AbstractMesh` types see [GeometryTypes.jl mesh documentation](http://juliageometry.github.io/GeometryTypes.jl/latest/types.html#Meshes-1).

### Meshing Algorithms

Three meshing algorithms exist:
* `MarchingCubes()`
* `MarchingTetrahedra()`
* `NaiveSurfaceNets()`

Each takes an optional `iso` and `eps` parameter, e.g. `MarchingCubes(0.0,1e-6)`.

Here `iso` controls the offset for the boundary detection. By default this is set to 0. `eps` is the detection tolerance for a voxel edge intersection.

Below is a comparison of the algorithms:

| Algorithm          | Accurate | Manifold | Performance Penalty | Face Type |
|--------------------|----------|----------|---------------------|-----------|
| MarchingCubes      | Yes      | No       | ~4x                 | Triangle  |
| MarchingTetrahedra | Yes      | Yes      | ~12x                | Triangle  |
| NaiveSurfaceNets   | No       | No       | 1x                  | Quad      |

## Credits

Marching Tetrahedra was originally implemented by @twadleigh in Meshes.jl.

## License

Available under the MIT "Expat" License, see [LICENSE.md](./LICENSE.md)
