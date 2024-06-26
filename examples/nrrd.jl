using Meshing
using FileIO
using NRRD
using WGLMakie
using Downloads
using GeometryBasics: GeometryBasics, Mesh

nrrd = Downloads.download("http://www.slicer.org/slicerWiki/images/0/00/CTA-cardio.nrrd")

# load the file as an AxisArray
ctacardio = load(nrrd)

# use marching cubes with isolevel at 100
#algo = MarchingCubes(iso=100)
# use marching tetrahedra with iso at 100
 algo = MarchingTetrahedra(iso=100)
# use Naive Surface Nets with iso at 100
# algo = NaiveSurfaceNets(iso=100)

# generate the mesh using marching cubes
vts, fcs = isosurface(ctacardio, algo)

WGLMakie.mesh(vts, map(v -> GeometryBasics.TriangleFace(v...), fcs))

