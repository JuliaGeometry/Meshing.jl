using Meshing
using Base.Test
using GeometryTypes


@testset "meshing" begin
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
            m1 = @test_nowarn marching_cubes(sdf, 0.0, GLNormalMesh)
            m2 = @test_nowarn GLNormalMesh(sdf, MarchingCubes())
            @test vertices(m1) == vertices(m2)
            @test faces(m1) == faces(m2)
        end

        @testset "marching tetrahedra" begin
            m1 = @test_warn "is deprecated" GLNormalMesh(sdf)
            m2 = @test_warn "is deprecated" GLNormalMesh(sdf, 1e-3)
            m3 = @test_nowarn GLNormalMesh(sdf, MarchingTetrahedra())
            @test vertices(m1) == vertices(m2) == vertices(m3)
            @test faces(m1) == faces(m2) == faces(m3)
        end
    end
end


if "--profile" in ARGS
    HomogenousMesh(s2)
    Profile.clear()
    @profile HomogenousMesh(s2)
    #ProfileView.view()
end
