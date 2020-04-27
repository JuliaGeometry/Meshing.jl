using BenchmarkTools
using Meshing
using GeometryBasics

# Define a parent BenchmarkGroup to contain our suite
const suite = BenchmarkGroup()
suite["SDF Construction"] = BenchmarkGroup()
suite["SDF Mesh"] = BenchmarkGroup()
suite["Function Mesh"] = BenchmarkGroup()

println("Benchmarking Meshing.jl...")

#
# Algorithms to benchmark
#

algos_sdf = [MarchingCubes(), MarchingTetrahedra(), NaiveSurfaceNets()]
algos_fn = [MarchingCubes(), MarchingTetrahedra(), NaiveSurfaceNets()]

#
# Benchmark SDF constructon for SDF and Direct sample comparisons
#

suite["SDF Construction"]["Torus"] = @benchmarkable SignedDistanceField(Rect(Vec(-2,-2,-2.),Vec(4,4,4.)),0.05) do v
    (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus
end

#
# Torus Constructions
#

sdf_torus = SignedDistanceField(Rect(Vec(-2,-2,-2.),Vec(4,4,4.)),0.05) do v
    (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus
end

fn_torus(v) = (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus

#
# Benchmark algorithms
#

for algo in algos_sdf
    suite["SDF Mesh"][string(typeof(algo))] = @benchmarkable HomogenousMesh(sdf_torus, $algo)
end

for algo in algos_fn
    suite["Function Mesh"][string(typeof(algo))] = @benchmarkable HomogenousMesh(fn_torus, Rect(Vec(-2,-2,-2.),Vec(4,4,4.)), (81,81,81), $algo)
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

#
# Perform benchmarks and print results
#

results = run(suite)

for trial in results
    ctx = IOContext(stdout, :verbose => true, :compact => false)
    println(ctx)
    println(ctx, trial.first)
    println(ctx, trial.second)
end
