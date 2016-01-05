using Meshing
using Base.Test
using GeometryTypes

# Produce a level set function that is a noisy version of the distance from
# the origin (such that level sets are noisy spheres).
#
# The noise should exercise marching tetrahedra's ability to produce a water-
# tight surface in all cases (unlike standard marching cubes).
#
N = 10
sigma = 1.0
distance = Float32[ sqrt(Float32(i*i+j*j+k*k)) for i = -N:N, j = -N:N, k = -N:N ]
distance = distance + sigma*rand(2*N+1,2*N+1,2*N+1)

# Extract an isosurface.
#
lambda = N-2*sigma # isovalue

msh = HomogenousMesh(distance,lambda)

s2 = SignedDistanceField(HyperRectangle(Vec(0,0,0.),Vec(1,1,1.))) do v
    sqrt(sum(dot(v,v))) - 1 # sphere
end

msh = HomogenousMesh(s2)
@test length(vertices(msh)) == 973
@test length(faces(msh)) == 1830

# Vertex interpolation
@test Meshing.vertex_interp(0, Vec(0,0,0), Vec(0,1,0), -1, 1) == Vec(0,0.5,0)
@test Meshing.vertex_interp(-1, Vec(0,0,0), Vec(0,1,0), -1, 1) == Vec(0,0,0)
@test Meshing.vertex_interp(1, Vec(0,0,0), Vec(0,1,0), -1, 1) == Vec(0,1,0)

# marching cubes
sdf = SignedDistanceField(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.))) do v
    sqrt(sum(dot(v,v))) - 1 # sphere
end

m = marching_cubes(sdf,0)
m2 = marching_cubes(sdf)
@test length(vertices(m)) == 10968
@test length(faces(m)) == 3656
@test m == m2

if "--profile" in ARGS
    HomogenousMesh(s2)
    Profile.clear()
    @profile HomogenousMesh(s2)
    #ProfileView.view()
end
