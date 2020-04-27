using Meshing
using GeometryBasics
using Juno
using Profile


# algo type to profile
algo = MarchingTetrahedra()

# set if using function or SDF variant
fn_mesh = false

## Constructions
torus = SignedDistanceField(Rect(Vec(-2,-2,-2.),Vec(4,4,4.)),0.05) do v
    (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus
end

torus_fn(v) = (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus

# clearout profile
Profile.clear()

# warmup and profile
if fn_mesh
    HomogenousMesh(torus_fn,Rect(Vec(-2,-2,-2.),Vec(4,4,4.)), algo) # use default 128 size for warmup
    Juno.@profiler HomogenousMesh(torus_fn,Rect(Vec(-2,-2,-2.),Vec(4,4,4.)), (512,512,512), algo)
else
    HomogenousMesh(torus, algo)
    Juno.@profiler HomogenousMesh(torus, algo)
end

Juno.profiler()
