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

"""
Generate a mesh using naive surface nets.
This takes the center of mass of the voxel as the vertex for each cube.
"""
function isosurface(sdf::AbstractArray{T, 3}, method::NaiveSurfaceNets, ::Type{VertType}=SVector{3,Float32}, ::Type{FaceType}=SVector{4, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2)) where {T, VertType, FaceType}

    scale = widths ./ VertType(size(sdf) .- 1)  # subtract 1 because an SDF with N points per side has N-1 cells

    dims = size(sdf)
    data = vec(sdf)

    # Run iso surface additions here
    # TODO
    if method.iso != 0.0
        for i = eachindex(data)
            data[i] -= method.iso
        end
    end

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

        yi=0
        @inbounds while yi<dims[2]-1

            xi=0
            @inbounds while xi < dims[1]-1

                # Read in 8 field values around this vertex and store them in an array
                # Also calculate 8-bit mask, like in marching cubes, so we can speed up sign checks later
                @inbounds grid = (data[n+1],
                                  data[n+2],
                                  data[n+dims[1]+1],
                                  data[n+dims[1]+2],
                                  data[n+dims[1]*2+1 + dims[1]*(dims[2]-2)],
                                  data[n+dims[1]*2+2 + dims[1]*(dims[2]-2)],
                                  data[n+dims[1]*3+1 + dims[1]*(dims[2]-2)],
                                  data[n+dims[1]*3+2 + dims[1]*(dims[2]-2)])

                mask = _get_cubeindex(grid, 0)

                # Check for early termination if cell does not intersect boundary
                if mask == 0x00 || mask == 0xff
                    xi += 1
                    n += 1
                    m += 1
                    continue
                end

                #Sum up edge intersections
                edge_mask = sn_edge_table[mask]

                _sn_add_verts!(xi, yi, zi, vertices, grid, edge_mask, buffer, m, scale, origin, method.eps, T, Val(true), VertType)

                #Now we need to add faces together, to do this we just loop over 3 basis components
                x = (xi,yi,zi)
                for i=0:2
                    #The first three entries of the edge_mask count the crossings along the edge
                    if (edge_mask & (1<<i)) == 0
                        continue
                    end

                    # i = axes we are point along.  iu, iv = orthogonal axes
                    iu = (i+1)%3
                    iv = (i+2)%3

                    #If we are on a boundary, skip it
                    if (x[iu+1] == 0 || x[iv+1] == 0)
                        continue
                    end

                    #Otherwise, look up adjacent edges in buffer
                    du = R[iu+1]
                    dv = R[iv+1]

                    #Remember to flip orientation depending on the sign of the corner.
                    if (mask & 0x01) != 0x00
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
                xi += 1
                n += 1
                m += 1
            end
            yi += 1
            n += 1
            m += 2
        end
        zi += 1
        n+=dims[1]
        buf_no = xor(buf_no,1)
        R[3]=-R[3]
    end
    #All done!  Return the result

    vertices, faces # faces are quads, indexed to vertices
end

"""
Generate a mesh using naive surface nets.
This takes the center of mass of the voxel as the vertex for each cube.
"""
function isosurface(f::Function, method::NaiveSurfaceNets,
                    ::Type{VertType}=SVector{3,Float32}, ::Type{FaceType}=SVector{4, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2),
                    samples::NTuple{3,T}=(50,50,50)) where {T <: Integer, VertType, FaceType}

    scale = widths ./ VertType(samples .- 1)  # subtract 1 because an SDF with N points per side has N-1 cells


    vertices = VertType[]
    faces = FaceType[]

    sizehint!(vertices,ceil(Int,maximum(samples)^2))
    sizehint!(faces,ceil(Int,maximum(samples)^2))

    n = 0
    R = Array{Int}([1, (samples[1]+1), (samples[1]+1)*(samples[2]+1)])
    buf_no = 1

    buffer = fill(zero(Int),R[3]*2)

    grid = Vector{eltype(VertType)}(undef,8)

    #March over the voxel grid
    zi = 0
    @inbounds while zi<samples[3]-1

        # m is the pointer into the buffer we are going to use.
        # This is slightly obtuse because javascript does not have good support for packed data structures, so we must use typed arrays :(
        # The contents of the buffer will be the indices of the vertices on the previous x/y slice of the volume
        m = 1 + (samples[1]+1) * (1 + buf_no * (samples[2]+1))

        yi=0
        @inbounds while yi<samples[2]-1

            xi=0
            @inbounds while xi < samples[1]-1

                # Read in 8 field values around this vertex and store them in an array
                points = (VertType(xi,yi,zi).* scale + origin,
                          VertType(xi+1,yi,zi).* scale + origin,
                          VertType(xi,yi+1,zi).* scale + origin,
                          VertType(xi+1,yi+1,zi).* scale + origin,
                          VertType(xi,yi,zi+1).* scale + origin,
                          VertType(xi+1,yi,zi+1).* scale + origin,
                          VertType(xi,yi+1,zi+1).* scale + origin,
                          VertType(xi+1,yi+1,zi+1).* scale + origin)

                if xi == 0
                    for i = 1:8
                        grid[i] = f(points[i])
                    end
                else
                    grid[1] = grid[2]
                    grid[2] = f(points[2])
                    grid[3] = grid[4]
                    grid[4] = f(points[4])
                    grid[5] = grid[6]
                    grid[6] = f(points[6])
                    grid[7] = grid[8]
                    grid[8] = f(points[8])
                end

                # Also calculate 8-bit mask, like in marching cubes, so we can speed up sign checks later
                mask = _get_cubeindex(grid, 0)

                # Check for early termination if cell does not intersect boundary
                if mask == 0x00 || mask == 0xff
                    xi += 1
                    n += 1
                    m += 1
                    continue
                end

                #Sum up edge intersections
                edge_mask = sn_edge_table[mask]

                # add vertices
                _sn_add_verts!(xi, yi, zi, vertices, grid, edge_mask, buffer, m, scale, origin, method.eps, eltype(VertType), Val(false), VertType)

                #Now we need to add faces together, to do this we just loop over 3 basis components
                x = (xi,yi,zi)
                for i=0:2
                    #The first three entries of the edge_mask count the crossings along the edge
                    if (edge_mask & (1<<i)) == 0
                        continue
                    end

                    # i = axes we are point along.  iu, iv = orthogonal axes
                    iu = (i+1)%3
                    iv = (i+2)%3

                    #If we are on a boundary, skip it
                    if (x[iu+1] == 0 || x[iv+1] == 0)
                        continue
                    end

                    #Otherwise, look up adjacent edges in buffer
                    du = R[iu+1]
                    dv = R[iv+1]

                    #Remember to flip orientation depending on the sign of the corner.
                    if (mask & 0x01) != 0x00
                        push!(faces,FaceType(buffer[m+1]+1, buffer[m-du+1]+1, buffer[m-du-dv+1]+1, buffer[m-dv+1]+1));
                    else
                        push!(faces,FaceType(buffer[m+1]+1, buffer[m-dv+1]+1, buffer[m-du-dv+1]+1, buffer[m-du+1]+1));
                    end
                end
                xi += 1
                n += 1
                m += 1
            end
            yi += 1
            n += 1
            m += 2
        end
        zi += 1
        n+=samples[1]
        buf_no = xor(buf_no,1)
        R[3]=-R[3]
    end
    #All done!  Return the result

    vertices, faces # faces are quads, indexed to vertices
end

@inline function _sn_add_verts!(xi, yi, zi, vertices, grid, edge_mask, buffer, m, scale, origin, eps, T, translate_pt, ::Type{VertType}) where {VertType}
    v = zero(VertType)
    e_count = 0

    #For every edge of the cube...
    @inbounds for i=0:11

        #Use edge mask to check if it is crossed
        if (edge_mask & (1<<i)) == 0
            continue
        end

        #If it did, increment number of edge crossings
        e_count += 1

        #Now find the point of intersection
        e0 = cube_edges[(i<<1)+1]       #Unpack vertices
        e1 = cube_edges[(i<<1)+2]
        g0 = grid[e0]                 #Unpack grid values
        g1 = grid[e1]
        t  = g0 - g1                 #Compute point of intersection
        if abs(t) > eps
            t = g0 / t
        else
            continue
        end

        #Interpolate vertices and add up intersections (this can be done without multiplying)
        # TODO lut table change may have made this incorrect
        # in which case e0=e0-1 and e1=e1-1
        xj, yj, zj = zero(T), zero(T), zero(T)
        a = e0 & 1
        b = e1 & 1
        (a != 0) && (xj += one(T))
        (a != b) && (xj += (a != 0 ? - t : t))
        a = e0 & 2
        b = e1 & 2
        (a != 0) && (yj += one(T))
        (a != b) && (yj += (a != 0 ? - t : t))
        a = e0 & 4
        b = e1 & 4
        (a != 0) && (zj += 1.0)
        (a != b) && (zj += (a != 0 ? - t : t))
        v += VertType(xj,yj,zj)

    end # edge check

    #Now we just average the edge intersections and add them to coordinate
    s = 1.0 / e_count
    if typeof(translate_pt) == Val{true}
        @inbounds v = (VertType(xi,yi,zi)  + s .* v) .* scale + origin
    else
        @inbounds v = (VertType(xi,yi,zi) + s .* v)# * scale[i] + origin[i]
    end

    #Add vertex to buffer, store pointer to vertex index in buffer
    buffer[m+1] = length(vertices)
    push!(vertices, v)
end

