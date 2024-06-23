# API

## Quick Start - isosurface

Given 3D levelset data such as a CT scan, we can do:

```julia
using Meshing

A = rand(50,50,50) # 3D Matrix

points,faces = isosurface(A)
```

An iso-level is specified within the algorithm specification as follows:

```julia
using Meshing

A = rand(50,50,50) # 3D Matrix

points,faces = isosurface(A, MarchingCubes(iso=1))
```

Alternatively, we can use the `isosurface` API to sample a function, avoiding allocations:

```julia
using Meshing
using LinearAlgebra
using StaticArrays

points, faces = isosurface(origin=SVector(-1,-1,-1.), widths = SVector(2,2,2.), samples = (40,40,40)) do v
    sqrt(sum(dot(v,v))) - 1
end

# by default MarchingCubes() is used, but we may specify a different algorithm as follows

points, faces = isosurface(MarchingTetrahedra(), origin=SVector(-1,-1,-1.), widths = SVector(2,2,2.), samples = (40,40,40)) do v
    sqrt(sum(dot(v,v))) - 1
end
```

## Isosurface

`isosurface` is the common and generic API for isosurface extraction with any type of abstract vector/vertex/face type.

```@docs
isosurface
```


## Meshing Algorithms

Three meshing algorithms exist:

* `MarchingCubes()`
* `MarchingTetrahedra()`
* `NaiveSurfaceNets()`

Each takes optional `iso`, `eps`, and `insidepositive` parameters, e.g. `MarchingCubes(iso=0.0,eps=1e-6,insidepositive=false)`.

Here `iso` controls the offset for the boundary detection. By default this is set to 0. `eps` is the detection tolerance for a voxel edge intersection.
`insidepositive` sets the sign convention for inside/outside the surface (default: false).

Users must construct an algorithm type and use it as an argument to a GeometryTypes mesh call or `isosurface` call.

Below is a comparison of the algorithms:

| Algorithm Type      | Face Type | Unique Vertices | Performance | Interpolation     |
|---------------------|-----------|-----------------|-------------|-------------------|
| Naive Surface Nets  | Quad      | Yes             | ~1x         | Voxel Edge Weight |
| Marching Cubes      | Triangle  | No/Partial      | 1x          | Linear on Edge    |
| Marching Tetrahedra | Triangle  | Yes             | 3x          | Linear on Edge    |

Visual Comparison:
From left: Marching Cubes, Naive Surface Nets, Marching Tetrahedra

![comparison](./img/comparison.png)

```@docs
MarchingCubes
MarchingTetrahedra
NaiveSurfaceNets
Meshing.AbstractMeshingAlgorithm
```
