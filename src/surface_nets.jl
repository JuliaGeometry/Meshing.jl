#
# SurfaceNets in Julia
#
# Ported from the Javascript implementation by Mikola Lysenko (C) 2012
# https://github.com/mikolalysenko/mikolalysenko.github.com/blob/master/Isosurface/js/surfacenets.js
# MIT License
#
# Based on: S.F. Gibson, "Constrained Elastic Surface Nets". (1998) MERL Tech Report.
#

#Look up Table
include("lut/sn.jl")

function isosurface(sdf::AbstractArray{T, 3}, method::NaiveSurfaceNets, ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{4, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2)) where {T, VertType, FaceType}

    scale = widths ./ VertType(size(sdf) .- 1)  # subtract 1 because an SDF with N points per side has N-1 cells

    dims = size(sdf)

    vertices = VertType[]
    faces = FaceType[]

    sizehint!(vertices,ceil(Int,maximum(dims)^2))
    sizehint!(faces,ceil(Int,maximum(dims)^2))

    n = 0
    R = Array{Int}([1, (dims[1]+1), (dims[1]+1)*(dims[2]+1)])
    buf_no = 1

    buffer = fill(zero(Int),R[3]*2)

    #March over the voxel grid
    zi = 0
    @inbounds while zi<dims[3]-1

        # m is the pointer into the buffer we are going to use.
        # This is slightly obtuse because javascript does not have good support for packed data structures, so we must use typed arrays :(
        # The contents of the buffer will be the indices of the vertices on the previous x/y slice of the volume
        m = 1 + (dims[1]+1) * (1 + buf_no * (dims[2]+1))

        for yi = 0:dims[2]-2
            for xi = 0:dims[1]-2

                inds = (xi,yi,zi)

                # Read in 8 field values around this vertex and store them in an array
                # Also calculate 8-bit mask, like in marching cubes, so we can speed up sign checks later
                @inbounds grid = (sdf[xi+1,yi+1,zi+1],
                                  sdf[xi+2,yi+1,zi+1],
                                  sdf[xi+1,yi+2,zi+1],
                                  sdf[xi+2,yi+2,zi+1],
                                  sdf[xi+1,yi+1,zi+2],
                                  sdf[xi+2,yi+1,zi+2],
                                  sdf[xi+1,yi+2,zi+2],
                                  sdf[xi+2,yi+2,zi+2])

                mask = method.insidepositive ? _get_cubeindex_pos(grid, method.iso) : _get_cubeindex(grid, method.iso)

                # Check for early termination if cell does not intersect boundary
                if _no_triangles(mask)
                    m += 1
                    continue
                end

                # iso level correction
                grid = grid .- method.iso

                #Sum up edge intersections
                edge_mask = sn_edge_table[mask]

                _sn_add_verts!(inds, vertices, grid, edge_mask, buffer, m, scale, origin, method.eps, VertType)

                #Now we need to add faces together, to do this we just loop over 3 basis components
                _sn_add_faces!(inds, faces, edge_mask, mask, buffer, m, R, FaceType)

                m += 1
            end
            m += 2
        end
        zi += 1
        buf_no = xor(buf_no,1)
        R[3] =- R[3]
    end

    return vertices, faces
end

function isosurface(f::Function, method::NaiveSurfaceNets,
                    ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{4, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2),
                    samples::NTuple{3,T}=_DEFAULT_SAMPLES) where {T <: Integer, VertType, FaceType}

    scale = widths ./ VertType(samples .- 1)  # subtract 1 because an SDF with N points per side has N-1 cells

    vertices = VertType[]
    faces = FaceType[]

    sizehint!(vertices,ceil(Int,maximum(samples)^2))
    sizehint!(faces,ceil(Int,maximum(samples)^2))

    R = Array{Int}([1, (samples[1]+1), (samples[1]+1)*(samples[2]+1)])
    buf_no = 1

    buffer = fill(zero(Int),R[3]*2)

    zv = zero(eltype(VertType))
    zvt = zero(VertType)
    grid = (zv,zv,zv,zv,zv,zv,zv,zv)
    points = (zvt,zvt,zvt,zvt,zvt,zvt,zvt,zvt)

    #March over the voxel grid
    zi = 0
    @inbounds while zi<samples[3]-1

        # m is the pointer into the buffer we are going to use.
        # This is slightly obtuse because javascript does not have good support for packed data structures, so we must use typed arrays :(
        # The contents of the buffer will be the indices of the vertices on the previous x/y slice of the volume
        m = 1 + (samples[1]+1) * (1 + buf_no * (samples[2]+1))

        for yi = 0:samples[2]-2
            for xi = 0:samples[1]-2

                inds = (xi,yi,zi)

                if xi == 0
                    points = (VertType(xi,yi,zi).* scale .+ origin,
                              VertType(xi+1,yi,zi).* scale .+ origin,
                              VertType(xi,yi+1,zi).* scale .+ origin,
                              VertType(xi+1,yi+1,zi).* scale .+ origin,
                              VertType(xi,yi,zi+1).* scale .+ origin,
                              VertType(xi+1,yi,zi+1).* scale .+ origin,
                              VertType(xi,yi+1,zi+1).* scale .+ origin,
                              VertType(xi+1,yi+1,zi+1).* scale .+ origin)
                    grid = (f(points[1]),
                            f(points[2]),
                            f(points[3]),
                            f(points[4]),
                            f(points[5]),
                            f(points[6]),
                            f(points[7]),
                            f(points[8]))
                else
                    points = (points[2],
                              VertType(xi+1,yi,zi).* scale .+ origin,
                              points[4],
                              VertType(xi+1,yi+1,zi).* scale .+ origin,
                              points[6],
                              VertType(xi+1,yi,zi+1).* scale .+ origin,
                              points[8],
                              VertType(xi+1,yi+1,zi+1).* scale .+ origin)
                    grid = (grid[2],
                            f(points[2]),
                            grid[4],
                            f(points[4]),
                            grid[6],
                            f(points[6]),
                            grid[8],
                            f(points[8]))
                end

                # Also calculate 8-bit mask, like in marching cubes, so we can speed up sign checks later
                mask = method.insidepositive ? _get_cubeindex_pos(grid, method.iso) : _get_cubeindex(grid, method.iso)

                # Check for early termination if cell does not intersect boundary
                if _no_triangles(mask)
                    m += 1
                    continue
                end

                # iso level correction
                grid = grid .- method.iso

                #Sum up edge intersections
                edge_mask = sn_edge_table[mask]

                # add vertices
                _sn_add_verts!(inds, vertices, grid, edge_mask, buffer, m, scale, origin, method.eps, VertType)

                #Now we need to add faces together, to do this we just loop over 3 basis components
                _sn_add_faces!(inds, faces, edge_mask, mask, buffer, m, R, FaceType)

                m += 1
            end
            m += 2
        end
        zi += 1
        buf_no = xor(buf_no,1)
        R[3]=-R[3]
    end

    return vertices, faces
end

@inline function _sn_add_verts!(inds, vertices, grid, edge_mask, buffer, m, scale, origin, eps, ::Type{VertType}) where {VertType}
    v = zero(VertType)
    T = eltype(VertType)
    e_count = 0

    #For every edge of the cube...
    @inbounds for i=0x00:0x0b

        #Use edge mask to check if it is crossed
        iszero(edge_mask & (0x0001<<i)) && continue

        #If it did, increment number of edge crossings
        e_count += 1

        #Now find the point of intersection
        e0 = cube_edges[(i<<0x01)+1]       #Unpack vertices
        e1 = cube_edges[(i<<0x01)+2]
        g0 = grid[e0+1]                 #Unpack grid values
        g1 = grid[e1+1]
        t  = g0 - g1                 #Compute point of intersection

        abs(t) <= eps && continue

        t = g0 / t

        #Interpolate vertices and add up intersections (this can be done without multiplying)
        a1, a2, a3 = 0x01 & e0, 0x02 & e0, 0x04 & e0
        b1, b2, b3 = 0x01 & e1, 0x02 & e1, 0x04 & e1
        xj = T(!iszero(a1))
        yj = T(!iszero(a2))
        zj = T(!iszero(a3))
        a1 != b1 && (xj += !iszero(a1) ? -t : t)
        a2 != b2 && (yj += !iszero(a2) ? -t : t)
        a3 != b3 && (zj += !iszero(a3) ? -t : t)
        v = v .+ VertType(xj,yj,zj)

    end # edge check

    #Now we just average the edge intersections and add them to coordinate
    s = one(T) / e_count
    v = (VertType(inds...)  .+ s .* v) .* scale .+ origin

    #Add vertex to buffer, store pointer to vertex index in buffer
    buffer[m+1] = length(vertices)
    push!(vertices, v)
end

function _sn_add_faces!(inds, faces, edge_mask, mask, buffer, m, R, ::Type{FaceType}) where {FaceType}
    for i = 0x00:0x02
        #The first three entries of the edge_mask count the crossings along the edge
        iszero(edge_mask & (0x0001<<i)) && continue

        # i = axes we are point along.  iu, iv = orthogonal axes
        iu = (i+0x01)%0x03
        iv = (i+0x02)%0x03

        #If we are on a boundary, skip it
        iszero(inds[iu+1]) || iszero(inds[iv+1]) && continue

        #Otherwise, look up adjacent edges in buffer
        du = R[iu+1]
        dv = R[iv+1]

        #Remember to flip orientation depending on the sign of the corner.
        if !iszero(mask & 0x01)
            if length(FaceType) == 4
                push!(faces,FaceType(buffer[m+1]+1, buffer[m-du+1]+1, buffer[m-du-dv+1]+1, buffer[m-dv+1]+1))
            elseif length(FaceType) == 3
                push!(faces,FaceType(buffer[m+1]+1, buffer[m-du+1]+1, buffer[m-du-dv+1]+1))
                push!(faces,FaceType(buffer[m-du-dv+1]+1, buffer[m-dv+1]+1, buffer[m+1]+1))
            end
        else
            if length(FaceType) == 4
                push!(faces,FaceType(buffer[m+1]+1, buffer[m-dv+1]+1, buffer[m-du-dv+1]+1, buffer[m-du+1]+1))
            elseif length(FaceType) == 3
                push!(faces,FaceType(buffer[m+1]+1, buffer[m-dv+1]+1, buffer[m-du-dv+1]+1))
                push!(faces,FaceType(buffer[m-du-dv+1]+1, buffer[m-du+1]+1, buffer[m+1]+1))
            end
        end
    end
end
