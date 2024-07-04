# Examples

## NRRD Data

The file for this example can be found here: [http://www.slicer.org/slicerWiki/images/0/00/CTA-cardio.nrrd](http://www.slicer.org/slicerWiki/images/0/00/CTA-cardio.nrrd)

```julia
using Meshing
using FileIO
using NRRD
using WGLMakie
using Downloads
using GeometryBasics

nrrd = Downloads.download("http://www.slicer.org/slicerWiki/images/0/00/CTA-cardio.nrrd")

# load the file as an AxisArray
ctacardio = load(nrrd)

# use marching cubes with isolevel at 100
#algo = MarchingCubes(iso=100)
# use marching tetrahedra with iso at 100
 algo = MarchingTetrahedra(iso=100)


# generate the mesh using marching cubes
vts, fcs = isosurface(ctacardio, algo)

WGLMakie.mesh(vts, map(v -> GeometryBasics.TriangleFace(v...), fcs))

```

![cta cardio](./img/ctacardio.png)

## Functions

```julia
using Meshing
using FileIO # MeshIO should also be installed
using GeometryBasics

gyroid(v) = cos(v[1])*sin(v[2])+cos(v[2])*sin(v[3])+cos(v[3])*sin(v[1])
gyroid_shell(v) = max(gyroid(v)-0.4,-gyroid(v)-0.4)

# generate directly using GeometryBasics API
# Rect specifies the sampling intervals
gy_mesh = Mesh(gyroid_shell, Rect(Vec(0,0,0),Vec(pi*4,pi*4,pi*4)),
                       MarchingCubes(), samples=(50,50,50))

save("gyroid.ply", gy_mesh)

# view with Makie
import Makie
using LinearAlgebra
Makie.mesh(gy_mesh, color=[norm(v) for v in coordinates(gy_mesh)])
```

![gyroid](./img/gyroid.png)


## Isocaps

We do not provide an equivalent to `isocaps` in Matlab, though
a similar result may be achieved by setting the boundary to a large value:

```julia
using GeometryTypes

gyroid(v) = cos(v[1])*sin(v[2])+cos(v[2])*sin(v[3])+cos(v[3])*sin(v[1])
gyroid_shell(v) = max(gyroid(v)-0.4,-gyroid(v)-0.4)

xr,yr,zr = ntuple(_->LinRange(0,pi*4,50),3)

A = [gyroid_shell((x,y,z)) for x in xr, y in yr, z in zr]
A[1,:,:] .= 1e10
A[:,1,:] .= 1e10
A[:,:,1] .= 1e10
A[end,:,:] .= 1e10
A[:,end,:] .= 1e10
A[:,:,end] .= 1e10

gy_mesh = GLNormalMesh(A, MarchingCubes())

# view with Makie
import Makie
using LinearAlgebra
Makie.mesh(gy_mesh, color=[norm(v) for v in gy_mesh.vertices])
```
