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
function surface_nets(data::Vector{T}, dims,eps,scale,origin) where {T}

    vertices = Point{3,T}[]
    faces = Face{4,Int}[]

    sizehint!(vertices,ceil(Int,maximum(dims)^2/2))
    sizehint!(faces,ceil(Int,maximum(dims)^2/2))

    n = 0
    R = Array{Int}([1, (dims[1]+1), (dims[1]+1)*(dims[2]+1)])
    buf_no = 1

    buffer = fill(zero(Int),R[3]*2)

    v = Vector{T}([0.0,0.0,0.0])

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
                mask = 0x00
                @inbounds grid = (data[n+1],
                                  data[n+2],
                                  data[n+dims[1]+1],
                                  data[n+dims[1]+2],
                                  data[n+dims[1]*2+1 + dims[1]*(dims[2]-2)],
                                  data[n+dims[1]*2+2 + dims[1]*(dims[2]-2)],
                                  data[n+dims[1]*3+1 + dims[1]*(dims[2]-2)],
                                  data[n+dims[1]*3+2 + dims[1]*(dims[2]-2)])

                signbit(grid[1]) && (mask |= 0x01)
                signbit(grid[2]) && (mask |= 0x02)
                signbit(grid[3]) && (mask |= 0x04)
                signbit(grid[4]) && (mask |= 0x08)
                signbit(grid[5]) && (mask |= 0x10)
                signbit(grid[6]) && (mask |= 0x20)
                signbit(grid[7]) && (mask |= 0x40)
                signbit(grid[8]) && (mask |= 0x80)

                # Check for early termination if cell does not intersect boundary
                if mask == 0x00 || mask == 0xff
                    xi += 1
                    n += 1
                    m += 1
                    continue
                end

                #Sum up edge intersections
                edge_mask = sn_edge_table[mask+1]
                v[1] = 0.0
                v[2] = 0.0
                v[3] = 0.0
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
                    g0 = grid[e0+1]                 #Unpack grid values
                    g1 = grid[e1+1]
                    t  = g0 - g1                 #Compute point of intersection
                    if abs(t) > eps
                      t = g0 / t
                    else
                      continue
                    end

                    #Interpolate vertices and add up intersections (this can be done without multiplying)
                    k = 1
                    for j = 1:3
                        a = e0 & k
                        b = e1 & k
                        (a != 0) && (v[j] += 1.0)
                        if a != b
                            v[j] += (a != 0 ? - t : t)
                        end
                        k<<=1
                    end
                end # edge check

                #Now we just average the edge intersections and add them to coordinate
                s = 1.0 / e_count
                @inbounds v[1] = (xi + s * v[1]) * scale[1] + origin[1]
                @inbounds v[2] = (yi + s * v[2]) * scale[2] + origin[2]
                @inbounds v[3] = (zi + s * v[3]) * scale[3] + origin[3]

                #Add vertex to buffer, store pointer to vertex index in buffer
                buffer[m+1] = length(vertices)
                push!(vertices, Point{3,T}(v[1],v[2],v[3]))

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
                        push!(faces,Face{4,Int}(buffer[m+1]+1, buffer[m-du+1]+1, buffer[m-du-dv+1]+1, buffer[m-dv+1]+1));
                    else
                        push!(faces,Face{4,Int}(buffer[m+1]+1, buffer[m-dv+1]+1, buffer[m-du-dv+1]+1, buffer[m-du+1]+1));
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
function surface_nets(f::Function, dims::NTuple{3,Int},eps,scale,origin)

    # TODO
    T = Float64

    vertices = Point{3,T}[]
    faces = Face{4,Int}[]

    sizehint!(vertices,ceil(Int,maximum(dims)^2/2))
    sizehint!(faces,ceil(Int,maximum(dims)^2/2))

    n = 0
    R = Array{Int}([1, (dims[1]+1), (dims[1]+1)*(dims[2]+1)])
    buf_no = 1

    buffer = fill(zero(Int),R[3]*2)

    v = Vector{T}([0.0,0.0,0.0])
    grid = Vector{Float64}(undef,8)

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
                points = (Point{3,Float64}(xi,yi,zi).* scale + origin,
                          Point{3,Float64}(xi+1,yi,zi).* scale + origin,
                          Point{3,Float64}(xi,yi+1,zi).* scale + origin,
                          Point{3,Float64}(xi+1,yi+1,zi).* scale + origin,
                          Point{3,Float64}(xi,yi,zi+1).* scale + origin,
                          Point{3,Float64}(xi+1,yi,zi+1).* scale + origin,
                          Point{3,Float64}(xi,yi+1,zi+1).* scale + origin,
                          Point{3,Float64}(xi+1,yi+1,zi+1).* scale + origin)

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
                mask = 0x00
                signbit(grid[1]) && (mask |= 0x01)
                signbit(grid[2]) && (mask |= 0x02)
                signbit(grid[3]) && (mask |= 0x04)
                signbit(grid[4]) && (mask |= 0x08)
                signbit(grid[5]) && (mask |= 0x10)
                signbit(grid[6]) && (mask |= 0x20)
                signbit(grid[7]) && (mask |= 0x40)
                signbit(grid[8]) && (mask |= 0x80)

                # Check for early termination if cell does not intersect boundary
                if mask == 0x00 || mask == 0xff
                    xi += 1
                    n += 1
                    m += 1
                    continue
                end

                #Sum up edge intersections
                edge_mask = sn_edge_table[mask+1]
                v[1] = 0.0
                v[2] = 0.0
                v[3] = 0.0
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
                    g0 = grid[e0+1]                 #Unpack grid values
                    g1 = grid[e1+1]
                    t  = g0 - g1                 #Compute point of intersection
                    if abs(t) > eps
                      t = g0 / t
                    else
                      continue
                    end

                    #Interpolate vertices and add up intersections (this can be done without multiplying)
                    k = 1
                    for j = 1:3
                        a = e0 & k
                        b = e1 & k
                        (a != 0) && (v[j] += 1.0)
                        if a != b
                            v[j] += (a != 0 ? - t : t)
                        end
                        k<<=1
                    end
                end # edge check

                #Now we just average the edge intersections and add them to coordinate
                s = 1.0 / e_count
                @inbounds v[1] = (xi + s * v[1])# * scale[i] + origin[i]
                @inbounds v[2] = (yi + s * v[2])# * scale[i] + origin[i]
                @inbounds v[3] = (zi + s * v[3])# * scale[i] + origin[i]

                #Add vertex to buffer, store pointer to vertex index in buffer
                buffer[m+1] = length(vertices)
                push!(vertices, Point{3,T}(v[1],v[2],v[3]))

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
                        push!(faces,Face{4,Int}(buffer[m+1]+1, buffer[m-du+1]+1, buffer[m-du-dv+1]+1, buffer[m-dv+1]+1));
                    else
                        push!(faces,Face{4,Int}(buffer[m+1]+1, buffer[m-dv+1]+1, buffer[m-du-dv+1]+1, buffer[m-du+1]+1));
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


struct NaiveSurfaceNets{T} <: AbstractMeshingAlgorithm
    iso::T
    eps::T
end

NaiveSurfaceNets(iso::T1=0.0, eps::T2=1e-3) where {T1, T2} = NaiveSurfaceNets{promote_type(T1, T2)}(iso, eps)

function (::Type{MT})(sdf::SignedDistanceField, method::NaiveSurfaceNets) where {MT <: AbstractMesh}
    bounds = sdf.bounds
    orig = origin(bounds)
    w = widths(bounds)
    scale = w ./ Point(size(sdf) .- 1)  # subtract 1 because an SDF with N points per side has N-1 cells

    d = vec(sdf.data)

    # Run iso surface additions here as to not
    # penalize surface net inner loops, and we are using a copy anyway
    if method.iso != 0.0
        for i = eachindex(d)
            d[i] -= method.iso
        end
    end

    vts, fcs = surface_nets(d,
                            size(sdf.data),
                            method.eps,
                            scale,
                            orig)
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function, bounds::HyperRectangle, size::NTuple{3,Int}, method::NaiveSurfaceNets) where {MT <: AbstractMesh}
    orig = origin(bounds)
    w = widths(bounds)
    scale = w ./ Point(size .- 1)  # subtract 1 because an SDF with N points per side has N-1 cells

    # TODO ISO val

    vts, fcs = surface_nets(f,
                            size,
                            method.eps,
                            scale,
                            orig)
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function, bounds::HyperRectangle, method::NaiveSurfaceNets;size::NTuple{3,Int}=(128,128,128)) where {MT <: AbstractMesh}
    MT(f,bounds,size,method)
end
