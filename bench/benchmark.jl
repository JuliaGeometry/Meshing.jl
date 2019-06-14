using BenchmarkTools
using Meshing
using GeometryTypes

# Define a parent BenchmarkGroup to contain our suite
const suite = BenchmarkGroup()

println("Benchmarking Meshing.jl...")

algos = [MarchingCubes(), MarchingTetrahedra(), NaiveSurfaceNets()]

torus = SignedDistanceField(HyperRectangle(Vec(-2,-2,-2.),Vec(4,4,4.)),0.05) do v
    (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus
end

for algo in algos
    suite[string(typeof(algo))] = @benchmarkable HomogenousMesh(torus, $algo)
end

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `suite` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals);
else
    tune!(suite)
    BenchmarkTools.save(paramspath, params(suite));
end

results = run(suite)
for trial in results
    ctx = IOContext(stdout, :verbose => true, :compact => false)
    println(ctx)
    println(ctx, trial.first)
    println(ctx, trial.second)
end
