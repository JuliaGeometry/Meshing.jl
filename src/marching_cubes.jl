
#Look up Table
include("lut/mc.jl")

function isosurface(sdf::AbstractArray{T,3}, method::MarchingCubes, X=-1:1, Y=-1:1, Z=-1:1) where {T}
    nx, ny, nz = size(sdf)

    vts = NTuple{3,float(T)}[]
    fcs = NTuple{3,Int}[]

    xp = LinRange(first(X), last(X), nx)
    yp = LinRange(first(Y), last(Y), ny)
    zp = LinRange(first(Z), last(Z), nz)

    @inbounds for xi = 1:nx-1, yi = 1:ny-1, zi = 1:nz-1

        iso_vals = (sdf[xi, yi, zi],
            sdf[xi+1, yi, zi],
            sdf[xi+1, yi+1, zi],
            sdf[xi, yi+1, zi],
            sdf[xi, yi, zi+1],
            sdf[xi+1, yi, zi+1],
            sdf[xi+1, yi+1, zi+1],
            sdf[xi, yi+1, zi+1])

        #Determine the index into the edge table which
        #tells us which vertices are inside of the surface
        cubeindex = _get_cubeindex(iso_vals, method.iso)

        # Cube is entirely in/out of the surface
        _no_triangles(cubeindex) && continue

        points = mc_vert_points(xi, yi, zi, xp, yp, zp)

        # process the voxel
        process_mc_voxel!(vts, fcs, cubeindex, points, method.iso, iso_vals)
    end
    vts, fcs
end


function isosurface(f::Function, method::MarchingCubes, X=-1:1, Y=-1:1, Z=-1:1; samples::NTuple{3,T}=_DEFAULT_SAMPLES) where {T<:Integer}

    nx, ny, nz = samples[1], samples[2], samples[3]

    # find widest type
    FT = promote_type(eltype(first(X)), eltype(first(Y)), eltype(first(Z)), eltype(T))

    vts = NTuple{3,float(FT)}[]
    fcs = NTuple{3,Int}[]

    xp = LinRange(first(X), last(X), nx)
    yp = LinRange(first(Y), last(Y), ny)
    zp = LinRange(first(Z), last(Z), nz)

    @inbounds for xi = 1:nx-1, yi = 1:ny-1, zi = 1:nz-1

        points = mc_vert_points(xi, yi, zi, xp, yp, zp)

        iso_vals = (f(points[1]),
            f(points[2]),
            f(points[3]),
            f(points[4]),
            f(points[5]),
            f(points[6]),
            f(points[7]),
            f(points[8]))

        #Determine the index into the edge table which
        #tells us which vertices are inside of the surface
        cubeindex = _get_cubeindex(iso_vals, method.iso)

        # Cube is entirely in/out of the surface
        _no_triangles(cubeindex) && continue

        process_mc_voxel!(vts, fcs, cubeindex, points, method.iso, iso_vals)
    end
    vts, fcs
end


function process_mc_voxel!(vts, fcs, cubeindex, points, iso, iso_vals)

    fct = length(vts)

    @inbounds begin
        # Add the vertices
        vert_to_add = _mc_verts[cubeindex]
        for i = 1:12
            vt = vert_to_add[i]
            iszero(vt) && break
            ed = _mc_edge_list[vt]
            push!(vts, vertex_interp(iso, points[ed[1]], points[ed[2]], iso_vals[ed[1]], iso_vals[ed[2]]))
        end

        # Add the faces
        offsets = _mc_connectivity[_mc_eq_mapping[cubeindex]]

        # There is atleast one face so we can push it immediately
        push!(fcs, (fct + 3, fct + 2, fct + 1))

        for i in (1, 4, 7, 10)
            iszero(offsets[i]) && return
            push!(fcs, (fct + offsets[i+2], fct + offsets[i+1], fct + offsets[i]))
        end
    end
end


"""    vertex_interp(iso, p1, p2, valp1, valp2)

Linearly interpolate the position where an isosurface cuts
an edge between two vertices, each with their own scalar value
"""
function vertex_interp(iso, p1, p2, valp1, valp2)
    mu = (iso - valp1) / (valp2 - valp1)
    p = p1 .+ mu .* (p2 .- p1)
    return p
end

"""
    mc_vert_points(xi,yi,zi,xp,yp,zp)

Returns a tuple of 8 points corresponding to each corner of a cube
"""
function mc_vert_points(xi, yi, zi, xp, yp, zp)
    ((xp[xi], yp[yi], zp[zi]),
        (xp[xi+1], yp[yi], zp[zi]),
        (xp[xi+1], yp[yi+1], zp[zi]),
        (xp[xi], yp[yi+1], zp[zi]),
        (xp[xi], yp[yi], zp[zi+1]),
        (xp[xi+1], yp[yi], zp[zi+1]),
        (xp[xi+1], yp[yi+1], zp[zi+1]),
        (xp[xi], yp[yi+1], zp[zi+1]))
end
