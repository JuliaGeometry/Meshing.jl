
using Meshing
using GeometryBasics
using LinearAlgebra: dot, norm
using FileIO

sphere(v) = sqrt(sum(v.^2)) - 1

sdf = SignedDistanceField(sphere, Rect(Vec(-1,-1,-1.), Vec(2,2,2.)))

mesh_type = HomogenousMesh{Point{3,Float64}, TriangleFace{Int}}

mc = mesh_type(sdf, MarchingCubes())
mt = mesh_type(sdf, MarchingTetrahedra())
ns = mesh_type(sdf, NaiveSurfaceNets())

mcf = mesh_type(sphere, sdf.bounds, MarchingCubes(), samples=size(sdf))
mtf = mesh_type(sphere, sdf.bounds, MarchingTetrahedra(), samples=size(sdf))
nsf = mesh_type(sphere, sdf.bounds, NaiveSurfaceNets(), samples=size(sdf))

# shift
@assert length(mc.vertices) == length(mcf.vertices)
@assert length(mt.vertices) == length(mtf.vertices)
@assert length(ns.vertices) == length(nsf.vertices)

for i in eachindex(mt.vertices)
    mt.vertices[i] = mt.vertices[i] .+ Point(3,0,0)
    mtf.vertices[i] = mtf.vertices[i] .+ Point(3,3,0)
end
for i in eachindex(mc.vertices)
    mc.vertices[i] = mc.vertices[i] .+ Point(-3,0,0)
    mcf.vertices[i] = mcf.vertices[i] .+ Point(-3,3,0)
end
for i in eachindex(ns.vertices)
    #ns.vertices[i] = ns.vertices[i] .+ Point(0,0,0)
    nsf.vertices[i] = nsf.vertices[i] .+ Point(0,3,0)
end

mesh = merge(mc,mt,ns,mcf,mtf,nsf)

# save all in one PLY file
save("output_comparison.ply",mesh)
