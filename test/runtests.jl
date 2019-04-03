using Meshing
using Test
using GeometryTypes
using ForwardDiff
using Profile
using Statistics: mean
using LinearAlgebra: dot, norm
using MeshIO
using FileIO
using BenchmarkTools


@testset "meshing" begin
    @testset "surface nets" begin
          sdf_sphere = SignedDistanceField(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.))) do v
              sqrt(sum(dot(v,v))) - 1 # sphere
          end
          sdf_torus = SignedDistanceField(HyperRectangle(Vec(-2,-2,-2.),Vec(4,4,4.)), 0.05) do v
              (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25
          end
          sphere = HomogenousMesh(sdf_sphere, NaiveSurfaceNets())
          torus = HomogenousMesh(sdf_torus, NaiveSurfaceNets())
          @test length(vertices(sphere)) == 1832
          @test length(vertices(torus)) == 5532
          @test length(faces(sphere)) == 1830
          @test length(faces(torus)) == 5532
          @test for vt in vertices(sphere)
              d = sqrt(sum(vt .^ 2))
              if 0.001 < (d-1) < 0.001
                  continue
              else
                  return false
              end
              true
          end
    end


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
        distance = distance + sigma*rand(2*N+1,2*N+1,2*N+1)

        # Extract an isosurface.
        #
        lambda = N-2*sigma # isovalue

        msh = HomogenousMesh(distance, MarchingTetrahedra(lambda))

        s2 = SignedDistanceField(HyperRectangle(Vec(0,0,0.),Vec(1,1,1.))) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        if "--profile" in ARGS
            HomogenousMesh(s2, MarchingTetrahedra())
            Profile.clear()
            @profile HomogenousMesh(s2, MarchingTetrahedra())
            #ProfileView.view()
        end

        msh = HomogenousMesh(s2, MarchingTetrahedra())
        @test length(vertices(msh)) == 973
        @test length(faces(msh)) == 1830
    end

    @testset "vertex interpolation" begin
        @test Meshing.vertex_interp(0, Vec(0,0,0), Vec(0,1,0), -1, 1) == Vec(0,0.5,0)
        @test Meshing.vertex_interp(-1, Vec(0,0,0), Vec(0,1,0), -1, 1) == Vec(0,0,0)
        @test Meshing.vertex_interp(1, Vec(0,0,0), Vec(0,1,0), -1, 1) == Vec(0,1,0)
    end

    @testset "marching cubes" begin
        sdf = SignedDistanceField(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.))) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        m = marching_cubes(sdf,0)
        m2 = marching_cubes(sdf)
        @test length(vertices(m)) == 10968
        @test length(faces(m)) == 3656
        @test m == m2

    end
    @testset "respect origin" begin
        # verify that when we construct a mesh, that mesh:
        #   a) respects the origin of the SDF
        #   b) is correctly spaced so that symmetric functions create symmetric meshes
        f = x -> norm(x)
        bounds = HyperRectangle(Vec(-1, -1, -1), Vec(2, 2, 2))
        resolution = 0.1
        sdf = SignedDistanceField(f, bounds, resolution)

        for algorithm in (MarchingCubes(0.5), MarchingTetrahedra(0.5))
            mesh = @inferred GLNormalMesh(sdf, algorithm)
            # should be centered on the origin
            @test mean(vertices(mesh)) ≈ [0, 0, 0] atol=0.15*resolution
            # and should be symmetric about the origin
            @test maximum(vertices(mesh)) ≈ [0.5, 0.5, 0.5]
            @test minimum(vertices(mesh)) ≈ [-0.5, -0.5, -0.5]
        end
    end

    @testset "AbstractMeshingAlgorithm interface" begin
        f = x -> norm(x) - 0.5
        bounds = HyperRectangle(Vec(-1, -1, -1), Vec(2, 2, 2))
        resolution = 0.1
        sdf = SignedDistanceField(f, bounds, resolution)

        @testset "marching cubes" begin
            @test_nowarn GLNormalMesh(sdf, MarchingCubes())
            @inferred GLNormalMesh(sdf, MarchingCubes())
        end

        @testset "marching tetrahedra" begin
            @test_nowarn GLNormalMesh(sdf, MarchingTetrahedra())
            @inferred GLNormalMesh(sdf, MarchingTetrahedra())
            @test_nowarn GLNormalMesh(sdf.data, MarchingTetrahedra(0.5))
            @inferred GLNormalMesh(sdf.data, MarchingTetrahedra(0.5))
        end
    end

    @testset "mixed types" begin
        # https://github.com/JuliaGeometry/Meshing.jl/issues/9
        samples = randn(10, 10, 10)
        m1 = HomogenousMesh(SignedDistanceField(
            HyperRectangle(Vec(0,0,0), Vec(1, 1, 1)),
            samples), MarchingTetrahedra())
        m2 = HomogenousMesh(SignedDistanceField(
            HyperRectangle(Vec(0,0,0), Vec(1, 1, 1)),
            Float32.(samples)), MarchingTetrahedra())
        @test all(isapprox.(vertices(m1), vertices(m2), rtol=1e-6))
        @test all(faces(m1) == faces(m2))

        @testset "forward diff" begin
            # Demonstrate that we can even take derivatives through the meshing algorithm
            f = x -> norm(x)
            bounds = HyperRectangle(Vec(-1, -1, -1), Vec(2, 2, 2))
            resolution = 0.1
            sdf = SignedDistanceField(f, bounds, resolution)

            function surface_distance_from_origin(isovalue)
                # Compute the mean distance of each vertex in the isosurface
                # mesh from the origin. This should return a value equal to the
                # isovalue and should have a derivative of 1.0 w.r.t. the
                # isovalue. This function is just meant to serve as an example
                # of some quantity you might want to differentiate in a mesh,
                # and has the benefit for testing of having a trivial expected
                # solution.
                mesh = HomogenousMesh(sdf, MarchingTetrahedra(isovalue))
                mean(norm.(vertices(mesh)))
            end

            isoval = 0.85
            @test surface_distance_from_origin(isoval) ≈ isoval atol=1e-2
            @test ForwardDiff.derivative(surface_distance_from_origin, isoval) ≈ 1 atol=1e-2
        end

        @testset "type stability" begin
            # verify that we don't lose type stability just by mixing arguments
            # of different types
            data = randn(5, 5, 5)
            iso = 0.2
            eps = 1e-4
            @inferred(Meshing.marchingTetrahedra(data, iso, eps, Int))
            @inferred(Meshing.marchingTetrahedra(Float32.(data), Float64(iso), Float16(eps), Int32))
            @inferred(Meshing.marchingTetrahedra(Float64.(data), Float32(iso), Float64(eps), Int64))
        end
    end
end
