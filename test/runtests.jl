using Meshing
using Test
using ForwardDiff
using Statistics: mean
using LinearAlgebra: dot, norm
using Random

const algos = (MarchingCubes, MarchingTetrahedra, NaiveSurfaceNets)

sphere_function(v) = sqrt(sum(dot(v, v))) - 1
torus_function(v) = (sqrt(v[1]^2 + v[2]^2) - 0.5)^2 + v[3]^2 - 0.25
@testset "meshing algorithms" begin
    for algo in (MarchingCubes, MarchingTetrahedra, NaiveSurfaceNets)
        @test algo(5) == algo{Float64}(5.0, 0.001, true, false)
        @test algo(0x00, 0x00) == algo{UInt8}(0x00, 0x00, true, false)
        @test algo{UInt16}(5, 0x00) == algo{UInt16}(0x0005, 0x0000, true, false)
        @test algo{Float32}(1) ==  algo{Float32}(1.0f0, 0.001f0, true, false)
        @test algo{Float32}(1.0, 0x00) == algo{Float32}(1.0f0, 0.0f0, true, false)
        @test algo{Float32}(iso=1, eps=0x00, insidepositive=true, reduceverts=false) ==
            algo{Float32}(1.0f0, 0.0f0, false, true)
    end
end

@testset "_get_cubeindex" begin
    iso_vals1 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    iso1 = 0.35
    @test Meshing._get_cubeindex(iso_vals1, iso1) == 0x07

    iso_vals2 = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
    iso2 = 0.75
    @test Meshing._get_cubeindex(iso_vals2, iso2) == 0x07

    iso_vals3 = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
    iso3 = 0.5
    @test Meshing._get_cubeindex(iso_vals3, iso3) == 0xe0
end

@testset "vertex interpolation" begin
    @test Meshing.vertex_interp(0, (0,0,0), (0,1,0), -1, 1) == (0,0.5,0)
    @test Meshing.vertex_interp(-1, (0,0,0), (0,1,0), -1, 1) == (0,0,0)
    @test Meshing.vertex_interp(1, (0,0,0), (0,1,0), -1, 1) == (0,1,0)
end

@testset "respect origin" begin
    # verify that when we construct a mesh, that mesh:
    # a) respects the origin of the SDF
    # b) is correctly spaced so that symmetric functions create symmetric meshes
    f = x -> norm(x)
    resolution = 0.1

    for algorithm in (MarchingCubes(0.5),
                      MarchingTetrahedra(0.5),
                      MarchingCubes(iso=0.5, reduceverts=false),
                      MarchingTetrahedra(iso=0.5, reduceverts=false))
        @testset "$(typeof(algorithm))" begin
            # Extract isosurface using a function
            points, faces = isosurface(f, algorithm, samples=(50, 50, 50))

            # should be centered on the origin
            @test mean(collect.(points)) ≈ [0, 0, 0] atol=0.15*resolution
            # and should be symmetric about the origin

            for i in 1:3
                @test maximum(map(p -> p[i], collect.(points))) ≈ 0.5 atol=0.001
                @test minimum(map(p -> p[i], collect.(points))) ≈ -0.5 atol=0.001
            end
        end
    end

    #TODO: SurfaceNets functional variant is bugged
    # Naive Surface Nets has no accuracy guarantee, and is a weighted sum
    # so a larger tolerance is needed for this one. In addition,
    # quad -> triangle conversion is not functioning correctly
    # see: https://github.com/JuliaGeometry/GeometryTypes.jl/issues/169
    #points, faces = isosurface(f, NaiveSurfaceNets(0.5), samples=(21, 21, 21))

    # should be centered on the origin
    #@test mean(points) ≈ [0, 0, 0] atol=0.15*resolution
    # and should be symmetric about the origin
    #@test maximum(points) ≈ [0.5, 0.5, 0.5] atol=0.2
    #@test minimum(points) ≈ [-0.5, -0.5, -0.5] atol=0.2
end
#@testset "surface nets" begin
#
#
#    """
#    Test the isosurface function with NaiveSurfaceNets and a function.
#
#    This code tests the `isosurface` function with the `NaiveSurfaceNets` algorithm and two different level set functions: `sphere_function` and `torus_function`. It checks that the number of vertices and faces in the resulting meshes match expected values.
#    """
#    # Test the isosurface function with NaiveSurfaceNets and a function
#    points_sphere, faces_sphere = isosurface(sphere_function, NaiveSurfaceNets(), origin=SVector(-1, -1, -1), widths=SVector(2, 2, 2), samples=(21, 21, 21))
#    points_torus, faces_torus = isosurface(torus_function, NaiveSurfaceNets(), origin=SVector(-2, -2, -2), widths=SVector(4, 4, 4), samples=(81, 81, 81))
#
#    # Add assertions to check the output
#    @test length(points_sphere) == 1832
#    @test length(faces_sphere) == 1830
#    @test length(points_torus) == 5532
#    @test length(faces_torus) == 5532
#end

@testset "noisy spheres" begin
    # Produce a level set function that is a noisy version of the distance from
    # the origin (such that level sets are noisy spheres).
    #
    # The noise should exercise marching tetrahedra's ability to produce a water-
    # tight surface in all cases (unlike standard marching cubes).
    #
    N = 10
    sigma = 1.0
    distance = Float32[ sqrt(Float32(i*i+j*j+k*k)) for i = -N:N, j = -N:N, k = -N:N ]
    distance = distance + sigma*rand(Random.MersenneTwister(0), 2*N+1,2*N+1,2*N+1)

    # Extract an isosurface.
    #
    lambda = N-2*sigma # isovalue

    points, faces = isosurface(distance, MarchingTetrahedra(lambda))

    @test length(points) == 3466
    @test length(faces) == 6928
end

@testset "marching cubes" begin
    algo = MarchingCubes()

    # Extract isosurfaces using MarchingCubes
    points_mf, faces_mf = isosurface(sphere_function, algo, samples=(21, 21, 21))

    # Test the number of vertices and faces
    @test length(points_mf) == 7320
    @test length(faces_mf) == 3656
end
@testset "marching tetrahedra" begin
    a1 = MarchingTetrahedra()
    a2 = MarchingTetrahedra(reduceverts=false)

    # Extract isosurfaces using MarchingTetrahedra
    points_mt1, faces_mt1 = isosurface(v -> sqrt(sum(dot(v, v))) - 1, a1, samples=(21, 21, 21))
    points_mt2, faces_mt2 = isosurface(v -> sqrt(sum(dot(v, v))) - 1, a2, samples=(21, 21, 21))

    # Test the number of vertices and faces
    @test length(points_mt1) == 5582
    @test length(points_mt2) == 33480
    @test length(faces_mt1) == 11160
    @test length(faces_mt2) == length(faces_mt1)

    # When zero set is not completely within the box
    points_mt1_partial, faces_mt1_partial = isosurface(v -> v[3] - 50, a1, 0:50, 0:50, 40:60, samples=(2, 2, 2))
    points_mt2_partial, faces_mt2_partial = isosurface(v -> v[3] - 50, a2, 0:50, 0:50, 40:60, samples=(2, 2, 2))

    # Test the number of vertices and faces
    @test length(points_mt1_partial) == 9
    @test length(points_mt2_partial) == 24
    @test length(faces_mt1_partial) == 8
    @test length(faces_mt2_partial) == length(faces_mt1_partial)
end

@testset "function/array" begin

    for algo in algos
        @testset "$algo" begin
            # Extract isosurface using a function
            points_fn, faces_fn = isosurface(torus_function, algo(), origin=SVector(-2, -2, -2), widths=SVector(4, 4, 4), samples=(81, 81, 81))

            # Extract isosurface using an array
            torus_array = [torus_function(SVector(x, y, z)) for x in -2:0.05:2, y in -2:0.05:2, z in -2:0.05:2]
            points_arr, faces_arr = isosurface(torus_array, algo(), origin=SVector(-2, -2, -2), widths=SVector(4, 4, 4))

            # Test that the vertices and faces are the same for both cases
            @test_broken all(points_fn .≈ points_arr)
            @test_broken faces_fn == faces_arr
        end
    end
end

@testset "mixed types" begin
    # https://github.com/JuliaGeometry/Meshing.jl/issues/9
    samples = randn(10, 10, 10)

    # Extract isosurfaces using MarchingTetrahedra
    points_mt1, faces_mt1 = isosurface(samples, MarchingTetrahedra(), origin=SVector(0, 0, 0), widths=SVector(1, 1, 1))
    points_mt2, faces_mt2 = isosurface(Float32.(samples), MarchingTetrahedra(), origin=SVector(0, 0, 0), widths=SVector(1, 1, 1))

    @test all(isapprox.(points_mt1, points_mt2, rtol=1e-6))
    @test all(faces_mt1 == faces_mt2)

    #=
    @testset "forward diff" begin
        # Demonstrate that we can even take derivatives through the meshing algorithm
        origin = SVector(-1, -1, -1)
        widths = SVector(2, 2, 2)
        resolution = 0.1

        for algo in (MarchingCubes, MarchingTetrahedra)
            @testset "$algo" begin
                let
                    function surface_distance_from_origin(isovalue)
                        # Compute the mean distance of each vertex in the isosurface
                        # mesh from the origin. This should return a value equal to the
                        # isovalue and should have a derivative of 1.0 w.r.t. the
                        # isovalue. This function is just meant to serve as an example
                        # of some quantity you might want to differentiate in a mesh,
                        # and has the benefit for testing of having a trivial expected
                        # solution.
                        points, _ = isosurface(norm, algo(isovalue,  reduce_vert=false), samples=(21, 21, 21))
                        mean(norm.(points))
                    end

                    isoval = 0.85
                    @test surface_distance_from_origin(isoval) ≈ isoval atol=1e-2
                    @test ForwardDiff.derivative(surface_distance_from_origin, isoval) ≈ 1 atol=1e-2
                end
            end
        end
    end

    =#

    @testset "type stability" begin
        # verify that we don't lose type stability just by mixing arguments
        # of different types
        data = randn(5, 5, 5)
        iso = 0.2
        eps = 1e-4
        for algo_ty in (MarchingTetrahedra, MarchingCubes,NaiveSurfaceNets)
            algo1 = algo_ty(iso, eps)
            algo2 = algo_ty(Float64(iso), Float16(eps))
            algo3 = algo_ty(Float32(iso), Float64(eps))
            @inferred(isosurface(data, algo1, SVector{3, Float16}, SVector{3,UInt32}))
            @inferred(isosurface(Float32.(data), algo2, SVector{3, Float32}, SVector{3,Int16}))
            @inferred(isosurface(Float64.(data), algo3, SVector{3, Float64}, SVector{3,UInt}))
        end
    end

    @testset "Float16" begin
        sdf_torus = [((sqrt(x^2 + y^2) - 0.5)^2 + z^2 - 0.25) for x in -2:0.05:2, y in -2:0.05:2, z in -2:0.05:2]
        points_nsn, faces_nsn = isosurface(Float16.(sdf_torus), NaiveSurfaceNets(Float16(0.0)), origin=SVector(Float16(-2), Float16(-2), Float16(-2)), widths=SVector(Float16(4), Float16(4), Float16(4)))
        points_mt, faces_mt = isosurface(Float16.(sdf_torus), MarchingTetrahedra(Float16(0.0)), origin=SVector(Float16(-2), Float16(-2), Float16(-2)), widths=SVector(Float16(4), Float16(4), Float16(4)))
        points_mc, faces_mc = isosurface(Float16.(sdf_torus), MarchingCubes(Float16(0.0)), origin=SVector(Float16(-2), Float16(-2), Float16(-2)), widths=SVector(Float16(4), Float16(4), Float16(4)))

        @test_broken typeof(points_nsn) == Vector{SVector{3, Float16}}
        @test typeof(faces_nsn) == Vector{SVector{4, Int}}
        @test_broken typeof(points_mt) == Vector{SVector{3, Float16}}
        @test typeof(faces_mt) == Vector{SVector{3, Int}}
        @test_broken typeof(points_mc) == Vector{SVector{3, Float16}}
        @test typeof(faces_mc) == Vector{SVector{3, Int}}
    end
    @testset "Integers" begin
        A = rand(Int, 10,10,10)
        for algo in algos
            @testset "$algo" begin
                p, t = isosurface(A, algo())
                @test typeof(p) <: Vector{SVector{3,Float64}}
                @test typeof(t) <: Vector{ SVector{ Meshing.default_face_length(algo()), Int } }
            end
        end
    end
    @testset "Defaults" begin
        A = rand(10,10,10)
        @test isosurface(A) == isosurface(A, MarchingCubes())
    end
end