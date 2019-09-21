
using Meshing
using GeometryTypes
using LinearAlgebra: dot, norm
using FileIO

f(v) = sqrt(sum(dot(v,v))) - 1
sdf = SignedDistanceField(f,HyperRectangle(Vec(-1,-1,-1.), Vec(2,2,2.)))

mc = HomogenousMesh(sdf, MarchingCubes())
mt = HomogenousMesh(sdf, MarchingTetrahedra())
ns = HomogenousMesh(sdf, NaiveSurfaceNets())

# save the Sphere as a PLY file
save("sphere_mc.ply",mc)
save("sphere_mt.ply",mt)
save("sphere_ns.ply",ns)
