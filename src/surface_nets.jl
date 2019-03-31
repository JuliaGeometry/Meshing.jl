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
const cube_edges = [ 0, 1, 0, 2, 0, 4, 1, 3, 1, 5, 2, 3, 2, 6, 3, 7, 4, 5, 4, 6, 5, 7, 6, 7 ]
const sn_edge_table = [
   0,
   7,
   25,
   30,
   98,
   101,
   123,
   124,
   168,
   175,
   177,
   182,
   202,
   205,
   211,
   212,
   772,
   771,
   797,
   794,
   870,
   865,
   895,
   888,
   940,
   939,
   949,
   946,
   974,
   969,
   983,
   976,
   1296,
   1303,
   1289,
   1294,
   1394,
   1397,
   1387,
   1388,
   1464,
   1471,
   1441,
   1446,
   1498,
   1501,
   1475,
   1476,
   1556,
   1555,
   1549,
   1546,
   1654,
   1649,
   1647,
   1640,
   1724,
   1723,
   1701,
   1698,
   1758,
   1753,
   1735,
   1728,
   2624,
   2631,
   2649,
   2654,
   2594,
   2597,
   2619,
   2620,
   2792,
   2799,
   2801,
   2806,
   2698,
   2701,
   2707,
   2708,
   2372,
   2371,
   2397,
   2394,
   2342,
   2337,
   2367,
   2360,
   2540,
   2539,
   2549,
   2546,
   2446,
   2441,
   2455,
   2448,
   3920,
   3927,
   3913,
   3918,
   3890,
   3893,
   3883,
   3884,
   4088,
   4095,
   4065,
   4070,
   3994,
   3997,
   3971,
   3972,
   3156,
   3155,
   3149,
   3146,
   3126,
   3121,
   3119,
   3112,
   3324,
   3323,
   3301,
   3298,
   3230,
   3225,
   3207,
   3200,
   3200,
   3207,
   3225,
   3230,
   3298,
   3301,
   3323,
   3324,
   3112,
   3119,
   3121,
   3126,
   3146,
   3149,
   3155,
   3156,
   3972,
   3971,
   3997,
   3994,
   4070,
   4065,
   4095,
   4088,
   3884,
   3883,
   3893,
   3890,
   3918,
   3913,
   3927,
   3920,
   2448,
   2455,
   2441,
   2446,
   2546,
   2549,
   2539,
   2540,
   2360,
   2367,
   2337,
   2342,
   2394,
   2397,
   2371,
   2372,
   2708,
   2707,
   2701,
   2698,
   2806,
   2801,
   2799,
   2792,
   2620,
   2619,
   2597,
   2594,
   2654,
   2649,
   2631,
   2624,
   1728,
   1735,
   1753,
   1758,
   1698,
   1701,
   1723,
   1724,
   1640,
   1647,
   1649,
   1654,
   1546,
   1549,
   1555,
   1556,
   1476,
   1475,
   1501,
   1498,
   1446,
   1441,
   1471,
   1464,
   1388,
   1387,
   1397,
   1394,
   1294,
   1289,
   1303,
   1296,
   976,
   983,
   969,
   974,
   946,
   949,
   939,
   940,
   888,
   895,
   865,
   870,
   794,
   797,
   771,
   772,
   212,
   211,
   205,
   202,
   182,
   177,
   175,
   168,
   124,
   123,
   101,
   98,
   30,
   25,
   7,
   0 ]

#Internal buffer, this may get resized at run time

function surface_nets(data, dims)
    buffer = Array{Int32}(undef,4096)

    vertices = []
    faces = []
    n = 0
    x = [0,0,0]
    R = Array{Int32}([1, (dims[1]+1), (dims[1]+1)*(dims[2]+1)])
    grid = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    buf_no = 1

    #Resize buffer if necessary
    if R[3]*2 > length(buffer)
        buffer = Array{Int32}(0,R[3]*2)
    end

    for i = 1:length(buffer)
      buffer[i] = 0
   end

    #March over the voxel grid
    x[3] = 0
    while x[3]<dims[3]-1

        # m is the pointer into the buffer we are going to use.
        # This is slightly obtuse because javascript does not have good support for packed data structures, so we must use typed arrays :(
        # The contents of the buffer will be the indices of the vertices on the previous x/y slice of the volume
        m = 1 + (dims[1]+1) * (1 + buf_no * (dims[2]+1))

        x[2]=0
        while x[2]<dims[2]-1

            x[1]=0
            while x[1] < dims[1]-1

                # Read in 8 field values around this vertex and store them in an array
                # Also calculate 8-bit mask, like in marching cubes, so we can speed up sign checks later
                mask = 0
                g = 0
                idx = n
                for k = 0:1
                    for j=0:1
                        for i=0:1
                            p = data[idx+1]
                            grid[g+1] = p
                            mask |= (p < 0) ? (1<<g) : 0
                            g += 1
                            idx += 1
                        end
                        idx += dims[1]-2
                    end
                    idx += dims[1]*(dims[2]-2)
                end

                # Check for early termination if cell does not intersect boundary
                if mask == 0 || mask == 0xff
                    x[1] += 1
                    continue
                end

                #Sum up edge intersections
                edge_mask = sn_edge_table[mask+1]
                v = [0.0,0.0,0.0]
                e_count = 0

                #For every edge of the cube...
                for i=0:11

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
                    if abs(t) > 1e-6
                      t = g0 / t
                    else
                      continue
                    end

                    #Interpolate vertices and add up intersections (this can be done without multiplying)
                    j = 0
                    k = 1
                    while j<3
                        a = e0 & k
                        b = e1 & k;
                        if a != b
                            v[j+1] += (a != 0 ? 1.0 - t : t)
                        else
                            v[j+1] += (a != 0 ? 1.0 : 0)
                        end
                        j+=1
                        k<<=1
                    end
                end # edge check

                #Now we just average the edge intersections and add them to coordinate
                s = 1.0 / e_count
                for i=1:3
                    v[i] = x[i] + s * v[i]
                end

                #Add vertex to buffer, store pointer to vertex index in buffer
                buffer[m+1] = length(vertices)
                push!(vertices, v);

                #Now we need to add faces together, to do this we just loop over 3 basis components
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
                    if (mask & 1) != 0
                        push!(faces,[buffer[m+1]+1, buffer[m-du+1]+1, buffer[m-du-dv+1]+1, buffer[m-dv+1]+1]);
                    else
                        push!(faces,[buffer[m+1]+1, buffer[m-dv+1]+1, buffer[m-du-dv+1]+1, buffer[m-du+1]+1]);
                    end
                end
                x[1] += 1
                n += 1
                m += 1
            end
            x[2] += 1
            n += 1
            m += 2
        end
        x[3] += 1
        n+=dims[1]
        buf_no = xor(buf_no,1)
        R[3]=-R[3]
    end
    #All done!  Return the result
    vts = Point{3,Float64}[]
    fcs = Face{4,Int}[]
    for face in faces
        push!(fcs, Face{4,Int}(face...))
    end
    for vt in vertices
        push!(vts, Point{3,Float64}(vt...))
    end
    @show vts, fcs
    return HomogenousMesh(vts, fcs)
end
