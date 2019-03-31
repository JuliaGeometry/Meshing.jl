#
# SurfaceNets in Julia
#
# Ported from the Javascript implementation by Mikola Lysenko (C) 2012
# https://github.com/mikolalysenko/mikolalysenko.github.com/blob/master/Isosurface/js/surfacenets.js
# MIT License
#
# Based on: S.F. Gibson, "Constrained Elastic Surface Nets". (1998) MERL Tech Report.
#



# Precompute edge table, like Paul Bourke does.
# This saves a bit of time when computing the centroid of each boundary cell
cube_edges = Array{Int32}(undef,24)
edge_table = Array{Int32}(undef,256)

# Initialize the cube_edges table
# This is just the vertex number of each cube
k = 0
for i=0:7
    for j in (1,2,4)
        p = i^j
        if i <= p
            cube_edges[k+1] = i
            k += 1
            cube_edges[k+1] = p
            k += 1
        end
    end
end

# Initialize the intersection table.
# This is a 2^(cube configuration) ->  2^(edge configuration) map
# There is one entry for each possible cube configuration, and the output is a 12-bit vector enumerating all edges crossing the 0-level.
for i=0:255
    em = 0;
    for j=0:2:22
        a = !!(i & (1<<cube_edges[j+1]))
        b = !!(i & (1<<cube_edges[j+2]))
        em = em | (a != b ? (1 << (j >> 1)) : 0)
    end
    edge_table[i+1] = em
end

#Internal buffer, this may get resized at run time
buffer = Array{Int32}(undef,4096)

function surface_nets(data, dims)

    vertices = []
    faces = []
    n = 0
    x = Array{Int32}(0, 3)
    R = Array{Int32}([1, (dims[1]+1), (dims[1]+1)*(dims[2]+1)])
    grid = Array{Float32}(0,8)
    buf_no = 1

    #Resize buffer if necessary
    if R[3]*2 > length(buffer)
        buffer = Array{Int32}(undef,R[3]*2)
    end

    #March over the voxel grid
    #for(x[2]=0; x[2]<dims[2]-1; ++x[2], n+=dims[0], buf_no ^= 1, R[2]=-R[2]) {
    while x[3]<dims[3]-1
        x[3] += 1
        n+=dims[1]
        buf_no ^= 1
        R[3]=-R[3]

        # m is the pointer into the buffer we are going to use.
        # This is slightly obtuse because javascript does not have good support for packed data structures, so we must use typed arrays :(
        # The contents of the buffer will be the indices of the vertices on the previous x/y slice of the volume
        m = 1 + (dims[1]+1) * (1 + buf_no * (dims[2]+1))

        x[2]=0
        while x[2]<dims[2]-1
            x[2] += 1
            n += 1
            m += 2
            x[1]=0
            while x[1] < dims[1]-1
                x[1] += 1
                n += 1
                m += 1

                # Read in 8 field values around this vertex and store them in an array
                # Also calculate 8-bit mask, like in marching cubes, so we can speed up sign checks later
                mask = 0
                g = 0
                idx = n
                for k = 0:1
                    idx += dims[1]*(dims[2]-2)
                    for j=0:1
                        idx += dims[1]-2
                        for i=0:1
                            g += 1
                            idx += 1
                            var p = data[idx]
                            grid[g] = p
                            mask |= (p < 0) ? (1<<g) : 0
                        end
                    end
                end

                # Check for early termination if cell does not intersect boundary
                if mask == 0 || mask == 0xff
                    continue
                end

                #Sum up edge intersections
                edge_mask = edge_table[mask]
                v = [0.0,0.0,0.0]
                e_count = 0

                #For every edge of the cube...
                for i=0:11

                    #Use edge mask to check if it is crossed
                    if !(edge_mask & (1<<i))
                        continue
                    end

                    #If it did, increment number of edge crossings
                    e_count += 1

                    #Now find the point of intersection
                    e0 = cube_edges[ i<<1 ]       #Unpack vertices
                    e1 = cube_edges[(i<<1)+1]
                    g0 = grid[e0]                 #Unpack grid values
                    g1 = grid[e1]
                    t  = g0 - g1                 #Compute point of intersection
                    if abs(t) > 1e-6
                      t = g0 / t
                    else
                      continue
                    end

                    #Interpolate vertices and add up intersections (this can be done without multiplying)
                    j = 0
                    k = 1
                    while j<3
                        j+=1
                        k<<=1
                        a = e0 & k
                        b = e1 & k;
                        if a != b
                            v[j+1] += a ? 1.0 - t : t
                        else
                            v[j+1] += a ? 1.0 : 0
                        end
                    end
                end # edge check

                #Now we just average the edge intersections and add them to coordinate
                s = 1.0 / e_count
                for i=0:2
                    v[i] = x[i] + s * v[i]
                end

                #Add vertex to buffer, store pointer to vertex index in buffer
                buffer[m] = length(vertices)
                push!(vertices, v);

                #Now we need to add faces together, to do this we just loop over 3 basis components
                for i=0:2
                    #The first three entries of the edge_mask count the crossings along the edge
                    if !(edge_mask & (1<<i))
                        continue
                    end

                    # i = axes we are point along.  iu, iv = orthogonal axes
                    iu = (i+1)%3
                    iv = (i+2)%3;

                    #If we are on a boundary, skip it
                    if (x[iu+1] == 0 || x[iv+1] == 0)
                        continue
                    end

                    #Otherwise, look up adjacent edges in buffer
                    du = R[iu+1]
                    dv = R[iv+1];

                    #Remember to flip orientation depending on the sign of the corner.
                    if (mask & 1)
                        push!(faces,[buffer[m+1], buffer[m-du+1], buffer[m-du-dv+1], buffer[m-dv+1]]);
                    else
                        push(faces,[buffer[m+1], buffer[m-dv+1], buffer[m-du-dv+1], buffer[m-du+1]]);
                    end
                end
            end
        end
    #All done!  Return the result
    return vertices, faces
end
