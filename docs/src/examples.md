# Examples


# NRRD Data

The file for this example can be found here: http://www.slicer.org/slicerWiki/images/0/00/CTA-cardio.nrrd

```

using FileIO
using NRRD
using Meshing
using MeshIO
using GeometryTypes

here = dirname(@__FILE__)

# load the file as an AxisArray
ctacardio = load(here*"/../data/CTA-cardio.nrrd")
q = -100 # isolevel


# flip sign on nrrd data since we need inside to be negative
# TODO
for i in eachindex(ctacardio.data)
    ctacardio.data[i] = -ctacardio.data[i]
end

# convert AxisArray to SignedDistanceField
ctasdf = SignedDistanceField(HyperRectangle(Vec(0,0,0), Vec(10,10,10)),ctacardio.data)

# generate the mesh using marching cubes
mc = HomogenousMesh{Point{3,Float32},Face{3,Int}}(ctasdf, MarchingCubes(q))

# save the file as a PLY file
save("ctacardio_mc.ply", mc)
```
