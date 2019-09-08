using Meshing
using Test
using GeometryTypes
using ForwardDiff
using StaticArrays
using Statistics: mean
using LinearAlgebra: dot, norm

isosdf(sdf::SignedDistanceField, algo, args...) = isosurface(sdf.data, algo, origin=origin(sdf.bounds), widths=widths(sdf.bounds), args...)

@testset "meshing" begin

    @testset "surface nets" begin
        sdf_sphere = SignedDistanceField(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.))) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end
        sdf_torus = SignedDistanceField(HyperRectangle(Vec(-2,-2,-2.),Vec(4,4,4.)), 0.05) do v
            (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25
        end

        snf_torus = isosurface(NaiveSurfaceNets(), samples=(81,81,81), origin=Vec(-2,-2,-2.), widths=Vec(4,4,4.)) do v
            (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25
        end

        snf_sphere = isosurface(NaiveSurfaceNets(), samples=(21,21,21), origin=Vec(-1,-1,-1.), widths=Vec(2,2,2.)) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        sphere_vts, sphere_fcs = isosurface(sdf_sphere.data, NaiveSurfaceNets(), origin=origin(sdf_sphere.bounds), widths=widths(sdf_sphere.bounds))
        torus_vts, torus_fcs = isosurface(sdf_torus.data, NaiveSurfaceNets(), origin=origin(sdf_torus.bounds), widths=widths(sdf_torus.bounds))
        @test length(sphere_vts) == 1832
        @test length(torus_vts) == 5532
        @test length(sphere_fcs) == 1830
        @test length(torus_fcs) == 5532

        @test length(sphere_vts) == length(snf_sphere[1])
        @test length(torus_vts) == length(snf_torus[1])
        @test length(sphere_fcs) == length(snf_sphere[2])
        @test length(torus_fcs) == length(snf_torus[2])
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

        msh = isosurface(distance, MarchingTetrahedra(lambda))

        s2 = SignedDistanceField(HyperRectangle(Vec(0,0,0.),Vec(1,1,1.))) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        msh_vts, msh_fcs = isosurface(s2.data, MarchingTetrahedra(), origin=origin(s2.bounds), widths=widths(s2.bounds))
        @test length(msh_vts) == 973
        @test length(msh_fcs) == 1830
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

        mf_vts, mf_fcs = isosurface(algo, samples=(21,21,21), origin=Vec(-1,-1,-1.), widths=Vec(2,2,2.)) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        mf_pos_vts, mf_pos_fcs = isosurface(algo_pos, samples=(21,21,21), origin=Vec(-1,-1,-1.), widths=Vec(2,2,2.)) do v
            -sqrt(sum(dot(v,v))) + 1 # sphere positive inside
        end

        mfrv_vts, mfrv_fcs = isosurface(MarchingCubes(reduceverts=false), samples=(21,21,21), origin=Vec(-1,-1,-1.), widths=Vec(2,2,2.)) do v
            sqrt(sum(dot(v,v))) - 1 # sphere
        end

        @test length(mfrv_vts) == 10968
        m_vts, m_fcs = isosurface(sdf.data,algo,origin=origin(sdf.bounds), widths=widths(sdf.bounds))
        m2_vts, m2_fcs = isosurface(sdf.data,algo,origin=origin(sdf.bounds), widths=widths(sdf.bounds))
        @test length(m_vts) == 7320
        @test length(m_fcs) == 3656
        @test length(mf_fcs) == length(mfrv_fcs)
        @test length(mf_fcs) == length(mf_pos_fcs)
        @test m_vts == m2_vts
        @test m_fcs == m2_fcs
        @test length(m_vts) == length(mf_vts)
        @test length(mf_vts) == length(mf_pos_vts)
        @test length(m_fcs) == length(mf_fcs)
    end

    @testset "respect origin" begin
        # verify that when we construct a mesh, that mesh:
        #   a) respects the origin of the SDF
        #   b) is correctly spaced so that symmetric functions create symmetric meshes
        f = x -> norm(x)
        bounds = HyperRectangle(Vec(-1, -1, -1.), Vec(2, 2, 2.))
        resolution = 0.1
        sdf = SignedDistanceField(f, bounds, resolution)

        for algorithm in (MarchingCubes(0.5), MarchingTetrahedra(0.5))
            vts, _ = @inferred isosurface(sdf.data, algorithm, Vec{3,Float32})
            # should be centered on the origin
            @test isapprox(mean(vts), SVector(0, 0, 0), atol=0.15*resolution)
            # and should be symmetric about the origin
            @test maximum(vts) ≈ [0.5, 0.5, 0.5]
            @test minimum(vts) ≈ [-0.5, -0.5, -0.5]
        end

        # Test functional Meshing
        vts, _ = @inferred isosurface(f, MarchingCubes(0.5), Point{3, Float64}, samples=(21,21,21))
        # should be centered on the origin
        @test isapprox(mean(vts), SVector(0, 0, 0), atol=0.15*resolution)
        # and should be symmetric about the origin
        @test maximum(vts) ≈ [0.5, 0.5, 0.5]
        @test minimum(vts) ≈ [-0.5, -0.5, -0.5]

        # Naive Surface Nets has no accuracy guarantee, and is a weighted sum
        # so a larger tolerance is needed for this one. In addition,
        # quad -> triangle conversion is not functioning correctly
        # see: https://github.com/JuliaGeometry/GeometryTypes.jl/issues/169
        vts, _ = @inferred isosurface(sdf.data, NaiveSurfaceNets(0.5), Point{3,Float64})
        # should be centered on the origin
        @test isapprox(mean(vts), SVector(0, 0, 0), atol=0.15*resolution)
        # and should be symmetric about the origin
        @test maximum(vts) ≈ [0.5, 0.5, 0.5] atol=0.2
        @test minimum(vts) ≈ [-0.5, -0.5, -0.5] atol=0.2
    end

    @testset "AbstractMeshingAlgorithm interface" begin
        f = x -> norm(x) - 0.5
        bounds = HyperRectangle(Vec(-1, -1, -1), Vec(2, 2, 2))
        resolution = 0.1
        sdf = SignedDistanceField(f, bounds, resolution)

        @testset "marching cubes" begin
            @test_nowarn isosurface(sdf.data, MarchingCubes())
            @inferred isosurface(sdf.data, MarchingCubes())
        end

        @testset "marching tetrahedra" begin
            @test_nowarn isosurface(sdf.data, MarchingTetrahedra())
            @inferred isosurface(sdf.data, MarchingTetrahedra())
            @test_nowarn isosurface(sdf.data, MarchingTetrahedra(0.5))
            @inferred isosurface(sdf.data, MarchingTetrahedra(0.5))
        end
        @testset "naive surface nets" begin
            @test_nowarn isosurface(sdf.data, NaiveSurfaceNets())
            @inferred isosurface(sdf.data, NaiveSurfaceNets())
        end
    end

    @testset "mixed types" begin
        # https://github.com/JuliaGeometry/Meshing.jl/issues/9
        samples = randn(10, 10, 10)
        m1_vts, m1_fcs = isosdf(SignedDistanceField(
            HyperRectangle(Vec(0,0,0), Vec(1, 1, 1)),
            samples), MarchingTetrahedra(), SVector{3,Float64})
        m2_vts, m2_fcs = isosdf(SignedDistanceField(
            HyperRectangle(Vec(0,0,0), Vec(1, 1, 1)),
            Float32.(samples)), MarchingTetrahedra(), SVector{3,Float32})
        @test all(isapprox.(m1_vts, m2_vts, rtol=1e-6))
        @test all(m1_fcs == m2_fcs)

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
                mesh_vts, mesh_fcs = isosurface(sdf.data, MarchingTetrahedra(isovalue), SVector{3,typeof(isovalue)})
                mean(norm.(mesh_vts))
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
            @inferred isosurface(data, MarchingTetrahedra(iso,eps), Point{3, Float16}, Face{3,UInt32})
            @inferred isosurface(Float32.(data), MarchingTetrahedra(iso,eps), Point{3, Float32}, Face{3,UInt16})
            @inferred isosurface(Float64.(data), MarchingTetrahedra(Float32(iso),Float64(eps)), Point{3, Float64}, Face{3,UInt})
        end
        @testset "Float16" begin
            sdf_torus = SignedDistanceField(HyperRectangle(Vec{3,Float16}(-2,-2,-2.),
                                                           Vec{3,Float16}(4,4,4.)),
                                            0.1, Float16) do v
                (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25
            end
            @test typeof(isosurface(sdf_torus.data,NaiveSurfaceNets(), SVector{3,Float16})) ==
                Tuple{Array{SVector{3,Float16},1},Array{SVector{4,Int},1}}
            m2 = isosdf(sdf_torus,MarchingTetrahedra())
            m3 = isosdf(sdf_torus,MarchingCubes())
        end
    end
end
