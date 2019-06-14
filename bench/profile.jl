using Meshing
using GeometryTypes
using Juno
using Profile


algo = NaiveSurfaceNets()

torus = SignedDistanceField(HyperRectangle(Vec(-2,-2,-2.),Vec(4,4,4.)),0.05) do v
    (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus
end

m2 = HomogenousMesh(torus, algo)

Profile.clear()
Juno.@profiler HomogenousMesh(torus, algo)

Juno.profiler()
