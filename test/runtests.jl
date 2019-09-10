using Meshing
using Test
using GeometryTypes
using ForwardDiff
using Profile
using Statistics: mean
using LinearAlgebra: dot, norm


@testset "meshing" begin
    @testset "type helpers" begin
        dt = Meshing._determine_types
        @test dt(HomogenousMesh) == (Point{3,Float64}, Face{3,Int64})
        @test dt(GLNormalMesh) == (Point{3,Float32}, Face{3,OffsetInteger{-1,UInt32}})
        @test dt(HomogenousMesh, Float16) == (Point{3,Float16}, Face{3,Int64})
        @test dt(HomogenousMesh{Point{3, Float64}, Face{3, UInt32}}) == (Point{3,Float64}, Face{3,UInt32})
        @test dt(HomogenousMesh{Point{3, Float32}, Face{3, UInt32}}, Float64) == (Point{3,Float32}, Face{3,UInt32})
        @test dt(HomogenousMesh, Float64, 4) == (Point{3,Float64}, Face{4,Int64})
    end

    @testset "surface nets" begin
        sdf_sphere = SignedDistanceField(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.))) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end
        sdf_torus = SignedDistanceField(HyperRectangle(Vec(-2,-2,-2.),Vec(4,4,4.)), 0.05) do v
            (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25
        end

        snf_torus = HomogenousMesh(HyperRectangle(Vec(-2,-2,-2.),Vec(4,4,4.)), (81,81,81), NaiveSurfaceNets()) do v
            (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25
        end

        snf_sphere = HomogenousMesh(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.)), (21,21,21), NaiveSurfaceNets()) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        # test convience constructors
        HomogenousMesh(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.)), NaiveSurfaceNets()) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end
        HomogenousMesh(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.)), NaiveSurfaceNets(), size=(5,5,5)) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        sphere = HomogenousMesh(sdf_sphere, NaiveSurfaceNets())
        torus = HomogenousMesh(sdf_torus, NaiveSurfaceNets())
        @test length(vertices(sphere)) == 1832
        @test length(vertices(torus)) == 5532
        @test length(faces(sphere)) == 1830
        @test length(faces(torus)) == 5532

        @test length(vertices(sphere)) == length(vertices(snf_sphere))
        @test length(vertices(torus)) == length(vertices(snf_torus))
        @test length(faces(sphere)) == length(faces(snf_sphere))
        @test length(faces(torus)) == length(faces(snf_torus))
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

        algo = MarchingCubes()
        algo_pos = MarchingCubes(insidepositive=true)

        sdf = SignedDistanceField(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.))) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        mf = SimpleMesh(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.)),algo, size=(21,21,21)) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        mf_pos = SimpleMesh(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.)),algo_pos, size=(21,21,21)) do v
            -sqrt(sum(dot(v,v))) + 1 # sphere positive inside
        end

        mfrv = SimpleMesh(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.)),MarchingCubes(reduceverts=false), size=(21,21,21)) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        # convience constructors
        SimpleMesh(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.)), MarchingCubes()) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end
        SimpleMesh(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.)), MarchingCubes(), size=(5,6,7)) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        @test length(vertices(mfrv)) == 10968
        m = SimpleMesh(sdf,algo)
        m2 = SimpleMesh(sdf,algo)
        @test length(vertices(m)) == 7320
        @test length(faces(m)) == 3656
        @test length(faces(mf)) == length(faces(mfrv))
        @test length(faces(mf)) == length(faces(mf_pos))
        @test m == m2
        @test length(vertices(m)) == length(vertices(mf))
        @test length(vertices(mf)) == length(vertices(mf_pos))
        @test length(faces(m)) == length(faces(mf))
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

        # Test functional Meshing
        mesh = @inferred GLNormalMesh(f, HyperRectangle(Vec(-1, -1, -1), Vec(2, 2, 2)), (21,21,21), MarchingCubes(0.5))
        # should be centered on the origin
        @test mean(vertices(mesh)) ≈ [0, 0, 0] atol=0.15*resolution
        # and should be symmetric about the origin
        @test maximum(vertices(mesh)) ≈ [0.5, 0.5, 0.5]
        @test minimum(vertices(mesh)) ≈ [-0.5, -0.5, -0.5]

        # Naive Surface Nets has no accuracy guarantee, and is a weighted sum
        # so a larger tolerance is needed for this one. In addition,
        # quad -> triangle conversion is not functioning correctly
        # see: https://github.com/JuliaGeometry/GeometryTypes.jl/issues/169
        mesh = @inferred GLNormalMesh(sdf, NaiveSurfaceNets(0.5))
        # should be centered on the origin
        @test mean(vertices(mesh)) ≈ [0, 0, 0] atol=0.15*resolution
        # and should be symmetric about the origin
        @test maximum(vertices(mesh)) ≈ [0.5, 0.5, 0.5] atol=0.2
        @test minimum(vertices(mesh)) ≈ [-0.5, -0.5, -0.5] atol=0.2
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
        @testset "naive surface nets" begin
            @test_nowarn GLNormalMesh(sdf, NaiveSurfaceNets())
            @inferred GLNormalMesh(sdf, NaiveSurfaceNets())
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
            for algo_ty in (MarchingTetrahedra, MarchingCubes)
                algo1 = algo_ty(iso, eps)
                algo2 = algo_ty(Float64(iso), Float16(eps))
                algo3 = algo_ty(Float32(iso), Float64(eps))
                @inferred(isosurface(data, algo1, Point{3, Float16}, Face{3,UInt32}))
                @inferred(isosurface(Float32.(data), algo2, Point{3, Float32}, Face{3,Int16}))
                @inferred(isosurface(Float64.(data), algo3, Point{3, Float64}, Face{3,UInt}))
            end
        end
        @testset "Float16" begin
            sdf_torus = SignedDistanceField(HyperRectangle(Vec{3,Float16}(-2,-2,-2.),
                                                           Vec{3,Float16}(4,4,4.)),
                                            0.1, Float16) do v
                (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25
            end
            @test typeof(HomogenousMesh(sdf_torus,NaiveSurfaceNets())) ==
                         PlainMesh{Float16,Face{4,Int}}
            m2 = HomogenousMesh(sdf_torus,MarchingTetrahedra())
            m3 = HomogenousMesh(sdf_torus,MarchingCubes())
        end
    end
end
