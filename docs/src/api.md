# API

## Quick Start

The easiest way to work with Meshing is with GeometryTypes.
This package extends the [mesh constructors](http://juliageometry.github.io/GeometryTypes.jl/latest/types.html#Meshes-1)
from [GeometryTypes.jl](https://github.com/JuliaGeometry/GeometryTypes.jl) for convience.

The algorithms operate on a `Function`, `AbstractArray`, or `SignedDistanceField` and output a concrete `AbstractMesh`. For example:

```
using Meshing
using GeometryTypes
using LinearAlgebra: dot, norm
using FileIO

# Mesh an equation of sphere in the Axis-Aligned Bounding box starting
# at -1,-1,-1 and widths of 2,2,2 using Marching Cubes
m = GLNormalMesh(HyperRectangle(Vec(-1,-1,-1.), Vec(2,2,2.)), MarchingCubes()) do v
    sqrt(sum(dot(v,v))) - 1
end

# save the Sphere as a PLY file
save("sphere.ply",m)
```

For a full listing of concrete `AbstractMesh` types see [GeometryTypes.jl mesh documentation](http://juliageometry.github.io/GeometryTypes.jl/latest/types.html#Meshes-1).

## Meshing Algorithms

Three meshing algorithms exist:
* `MarchingCubes()`
* `MarchingTetrahedra()`
* `NaiveSurfaceNets()`

Each takes an optional `iso` and `eps` parameter, e.g. `MarchingCubes(0.0,1e-6)`.

Here `iso` controls the offset for the boundary detection. By default this is set to 0. `eps` is the detection tolerance for a voxel edge intersection.

Users must construct an algorithm type and use it as an argument to a GeometryTypes mesh call or `isosurface` call.

Below is a comparison of the algorithms:

|                     | Face Type | Unique Vertices | Performance | Interpolation     |
|---------------------|-----------|-----------------|-------------|-------------------|
| Naive Surface Nets  | Quad      | Yes             | ~1x         | Voxel Edge Weight |
| Marching Cubes      | Triangle  | No/Partial      | 1x          | Linear on Edge    |
| Marching Tetrahedra | Triangle  | Yes             | 3x          | Linear on Edge    |

```@docs
MarchingCubes
MarchingTetrahedra
NaiveSurfaceNets
Meshing.AbstractMeshingAlgorithm
```

## Isosurface

`isosurface` is the common and generic API for isosurface extraction with any type of abstract vector/vertex/face type.

```@docs
isosurface
```

## GeometryTypes

Meshing extends the mesh types in GeometryTypes for convience and use with visualization tools such as Makie and MeshCat.
