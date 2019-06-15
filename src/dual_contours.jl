#Cardinal directions
const dirs = ( Point(1,0,0), Point(0,1,0), Point(0,0,1) )

#Vertices of cube
const cube_verts = (Point(0, 0, 0), Point(0, 0, 1), Point(0, 1, 0),
                    Point(0, 1, 1), Point(1, 0, 0), Point(1, 0, 1),
                    Point(1, 1, 0), Point(1, 1, 1))


#Edges of cube
const cube_edges_dc = ((0, 1), (0, 2), (0, 1), (0, 4), (0, 2), (0, 4), (2, 3), (1, 3),
              (4, 5), (1, 5), (4, 6), (2, 6), (4, 5), (4, 6), (2, 3), (2, 6),
              (1, 3), (1, 5), (6, 7), (5, 7), (6, 7), (3, 7), (5, 7), (3, 7))


#Use non-linear root finding to compute intersection point
function estimate_hermite(f, v0, v1)
    l(t) = f((1.0-t)*v0 + t*v1)
    dl(t) = ForwardDiff.derivative(l,t)
    t0 = find_zero((l,dl),(0, 1))
    x0 = (1.0-t0)*v0 + t0*v1
    return (x0, ForwardDiff.gradient(f,x0))
end

#Input:
# f = implicit function
# df = gradient of f
# nc = resolution
function dual_contours(f, bounds::HyperRectangle, nc::NTuple{3,Int})

    orig = origin(bounds)
    width = widths(bounds)
    scale = width ./ Point(nc)
    #Compute vertices
    dc_verts = []
    vindex   = Dict()
    for x in 0:nc[1], y in 0:nc[2], z in 0:nc[3]
        idx = Point(x,y,z)
        o = Point(x,y,z) .* scale + orig

        #Get signs for cube
        cube_signs = [ f(o+(v.*scale))>0 for v in cube_verts ]

        if all(cube_signs) || !any(cube_signs)
            continue
        end

        #Estimate hermite data
        h_data = [ estimate_hermite(f, o+cube_verts[e[1]+1].*scale, o+cube_verts[e[2]+1].*scale)
            for e in cube_edges_dc if cube_signs[e[1]+1] != cube_signs[e[2]+1] ]

        #Solve qef to get vertex
        A = Array{Float64,2}(undef,length(h_data),3)
        for i in eachindex(h_data)
            A[i,:] = h_data[i][2]
        end
        b = [ dot(d[1],d[2]) for d in h_data ]
        v = A\b

        #Throw out failed solutions
        if norm(v-o) > 2
            continue
        end

        #Emit one vertex per every cube that crosses
        push!(vindex, idx => length(dc_verts))
        push!(dc_verts, (v, ForwardDiff.gradient(f,v)))
    end

    #Construct faces
    dc_faces = Face[]
    for x in 0:nc[1], y in 0:nc[2], z in 0:nc[3]

        idx = Point(x,y,z)
        if !haskey(vindex,idx)
            continue
        end

        #Emit one face per each edge that crosses
        o = Point(x,y,z) .* scale + orig
        for i in (1,2,3)
            for j in 1:i
                if haskey(vindex,idx+dirs[i]) && haskey(vindex,idx + dirs[j]) && haskey(vindex,idx + dirs[i] + dirs[j])
                    # determine orientation of the face from the true normal
                    v1, tn1 = dc_verts[vindex[idx]+1]
                    v2, tn2 = dc_verts[vindex[idx+dirs[i]]+1]
                    v3, tn3 = dc_verts[vindex[idx+dirs[j]]+1]

                    e1 = v1-v2
                    e2 = v1-v3
                    c = cross(e1,e2)
                    if dot(c, tn1) > 0
                        push!(dc_faces, Face{3,Int}(vindex[idx]+1, vindex[idx+dirs[i]]+1, vindex[idx+dirs[j]]+1) )
                        push!(dc_faces, Face{3,Int}(vindex[idx+dirs[i]+dirs[j]]+1, vindex[idx+dirs[j]]+1, vindex[idx+dirs[i]]+1) )
                    else
                        push!(dc_faces, Face{3,Int}(vindex[idx]+1, vindex[idx+dirs[j]]+1, vindex[idx+dirs[i]]+1) )
                        push!(dc_faces, Face{3,Int}(vindex[idx+dirs[i]+dirs[j]]+1, vindex[idx+dirs[i]]+1, vindex[idx+dirs[j]]+1) )
                    end
                end
            end
        end

    end
    return HomogenousMesh([Point(v[1]...) for v in dc_verts], dc_faces)
end

