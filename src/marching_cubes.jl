
#Look up Table
include("lut/mc.jl")

function isosurface(sdf::AbstractArray{T, 3}, method::MarchingCubes, ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{3, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2)) where {T, VertType, FaceType}
    nx, ny, nz = size(sdf)

    # we subtract one from the length along each axis because
    # an NxNxN SDF has N-1 cells on each axis
    s = VertType(widths[1]/(nx-1), widths[2]/(ny-1), widths[3]/(nz-1))

    # arrays for vertices and faces
    vts = VertType[]
    fcs = FaceType[]
    mt = max(nx,ny,nz)
    method.reduceverts && sizehint!(vts, mt*mt*5)
    !method.reduceverts && sizehint!(vts, mt*mt*6)
    sizehint!(fcs, mt*mt*2)
    @inbounds for zi = 1:nz-1, yi = 1:ny-1, xi = 1:nx-1


        iso_vals = (sdf[xi,yi,zi],
                    sdf[xi+1,yi,zi],
                    sdf[xi+1,yi+1,zi],
                    sdf[xi,yi+1,zi],
                    sdf[xi,yi,zi+1],
                    sdf[xi+1,yi,zi+1],
                    sdf[xi+1,yi+1,zi+1],
                    sdf[xi,yi+1,zi+1])

        #Determine the index into the edge table which
        #tells us which vertices are inside of the surface
        cubeindex = method.insidepositive ? _get_cubeindex_pos(iso_vals, method.iso) : _get_cubeindex(iso_vals, method.iso)

        # Cube is entirely in/out of the surface
        _no_triangles(cubeindex) && continue

        points = mc_vert_points(xi,yi,zi,s,origin,VertType)

        # Create the triangle
        method.reduceverts && _mc_unique_triangles!(method, points, iso_vals, vts, fcs, cubeindex, FaceType)
        !method.reduceverts && _mc_create_triangles!(method, points, iso_vals, vts, fcs, cubeindex, FaceType)
    end
    vts,fcs
end


function isosurface(f::Function, method::MarchingCubes, ::Type{VertType}=SVector{3,Float64}, ::Type{FaceType}=SVector{3, Int};
                    origin=VertType(-1,-1,-1), widths=VertType(2,2,2), samples::NTuple{3,T}=_DEFAULT_SAMPLES) where {T <: Integer, VertType, FaceType}

    nx, ny, nz = samples[1], samples[2], samples[3]

    # we subtract one from the length along each axis because
    # an NxNxN SDF has N-1 cells on each axis
    s = VertType(widths[1]/(nx-1), widths[2]/(ny-1), widths[3]/(nz-1))

    # arrays for vertices and faces
    vts = VertType[]
    fcs = FaceType[]
    mt = max(nx,ny,nz)
    method.reduceverts && sizehint!(vts, mt*mt*5)
    !method.reduceverts && sizehint!(vts, mt*mt*6)
    sizehint!(fcs, mt*mt*2)
    zv = zero(eltype(VertType))
    iso_vals = (zv,zv,zv,zv,zv,zv,zv,zv)
    @inbounds for zi = 1:nz-1, yi = 1:ny-1, xi = 1:nx-1

        points = mc_vert_points(xi,yi,zi,s,origin,VertType)

        if xi == 1
            iso_vals = (f(points[1]),
                        f(points[2]),
                        f(points[3]),
                        f(points[4]),
                        f(points[5]),
                        f(points[6]),
                        f(points[7]),
                        f(points[8]))
        else
            iso_vals = (iso_vals[2],
                        f(points[2]),
                        f(points[3]),
                        iso_vals[3],
                        iso_vals[6],
                        f(points[6]),
                        f(points[7]),
                        iso_vals[7])
        end

        #Determine the index into the edge table which
        #tells us which vertices are inside of the surface
        cubeindex = method.insidepositive ? _get_cubeindex_pos(iso_vals, method.iso) : _get_cubeindex(iso_vals, method.iso)

        # Cube is entirely in/out of the surface
        _no_triangles(cubeindex) && continue

        # Create the triangle
        method.reduceverts && _mc_unique_triangles!(method, points, iso_vals, vts, fcs, cubeindex, FaceType)
        !method.reduceverts && _mc_create_triangles!(method, points, iso_vals, vts, fcs, cubeindex, FaceType)

    end
    vts,fcs
end


"""
    _mc_create_triangles!(vts, fcs, vertlist, cubeindex, FaceType)

Create triangles by adding every point within a triangle to the vertex vector.
"""
@inline function _mc_create_triangles!(method, points, iso_vals, vts, fcs, cubeindex, ::Type{FaceType}) where {FaceType}

    vertlist = find_vertices_interp(points, iso_vals, cubeindex, method.iso, method.eps)

    fct = length(vts) + 3

    push!(vts, vertlist[tri_table[cubeindex][1]],
               vertlist[tri_table[cubeindex][2]],
               vertlist[tri_table[cubeindex][3]])
    push!(fcs, FaceType(fct, fct-1, fct-2))

    iszero(tri_table[cubeindex][4]) && return
    fct += 3
    push!(vts, vertlist[tri_table[cubeindex][4]],
               vertlist[tri_table[cubeindex][5]],
               vertlist[tri_table[cubeindex][6]])
    push!(fcs, FaceType(fct, fct-1, fct-2))

    iszero(tri_table[cubeindex][7]) && return
    fct += 3
    push!(vts, vertlist[tri_table[cubeindex][7]],
               vertlist[tri_table[cubeindex][8]],
               vertlist[tri_table[cubeindex][9]])
    push!(fcs, FaceType(fct, fct-1, fct-2))

    iszero(tri_table[cubeindex][10]) && return
    fct += 3
    push!(vts, vertlist[tri_table[cubeindex][10]],
               vertlist[tri_table[cubeindex][11]],
               vertlist[tri_table[cubeindex][12]])
    push!(fcs, FaceType(fct, fct-1, fct-2))

    iszero(tri_table[cubeindex][13]) && return
    fct += 3
    push!(vts, vertlist[tri_table[cubeindex][13]],
               vertlist[tri_table[cubeindex][14]],
               vertlist[tri_table[cubeindex][15]])
    push!(fcs, FaceType(fct, fct-1, fct-2))
end

"""
    _mc_unique_triangles!(vts, fcs, vertlist, cubeindex, FaceType)

Create triangles by only adding unique vertices within the voxel.
Each face may share a reference to a vertex with another face.
"""
function _mc_unique_triangles!(method, points, iso_vals, vts, fcs, cubeindex, ::Type{FaceType}) where {FaceType}
    @inbounds begin
        fct = length(vts)

        find_vertices_interp!(vts, points, iso_vals, cubeindex, method.iso, method.eps)

        offsets = _mc_connectivity[_mc_eq_mapping[cubeindex]]

        # There is atleast one face so we can push it immediately
        push!(fcs, FaceType(fct+offsets[3], fct+offsets[2], fct+offsets[1]))

        for i in (4,7,10,13)
            iszero(offsets[i]) && return
            push!(fcs, FaceType(fct+offsets[i+2], fct+offsets[i+1], fct+offsets[i]))
        end
    end
end

"""
    find_vertices_interp(points, iso_vals, cubeindex, iso, eps)

Find the vertices where the surface intersects the cube
"""
@inline function find_vertices_interp(points, iso_vals, cubeindex, iso, eps)
    VT = eltype(points)
    zv = Ref{VT}()[]
    (!iszero(edge_table[cubeindex] & 0x001) ? vertex_interp(iso,points[1],points[2],iso_vals[1],iso_vals[2], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x002) ? vertex_interp(iso,points[2],points[3],iso_vals[2],iso_vals[3], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x004) ? vertex_interp(iso,points[3],points[4],iso_vals[3],iso_vals[4], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x008) ? vertex_interp(iso,points[4],points[1],iso_vals[4],iso_vals[1], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x010) ? vertex_interp(iso,points[5],points[6],iso_vals[5],iso_vals[6], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x020) ? vertex_interp(iso,points[6],points[7],iso_vals[6],iso_vals[7], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x040) ? vertex_interp(iso,points[7],points[8],iso_vals[7],iso_vals[8], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x080) ? vertex_interp(iso,points[8],points[5],iso_vals[8],iso_vals[5], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x100) ? vertex_interp(iso,points[1],points[5],iso_vals[1],iso_vals[5], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x200) ? vertex_interp(iso,points[2],points[6],iso_vals[2],iso_vals[6], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x400) ? vertex_interp(iso,points[3],points[7],iso_vals[3],iso_vals[7], eps) : zv,
     !iszero(edge_table[cubeindex] & 0x800) ? vertex_interp(iso,points[4],points[8],iso_vals[4],iso_vals[8], eps) : zv)
end

"""
    find_vertices_interp(points, iso_vals, cubeindex, iso, eps)

Find the vertices where the surface intersects the cube
"""
function find_vertices_interp!(vts, points, iso_vals, cubeindex, iso, eps)
    @inbounds begin
        vert_to_add = _mc_verts[cubeindex]
        for i = 1:12
            vt = vert_to_add[i]
            iszero(vt) && break
            ed = _mc_edge_list[vt]
            push!(vts, vertex_interp(iso,points[ed[1]],points[ed[2]],iso_vals[ed[1]],iso_vals[ed[2]], eps))
        end
    end
end

"""
    vertex_interp(iso, p1, p2, valp1, valp2, eps = 0.00001)

Linearly interpolate the position where an isosurface cuts
an edge between two vertices, each with their own scalar value
"""
function vertex_interp(iso, p1, p2, valp1, valp2, eps = 0.00001)

    abs(iso - valp1) < eps && return p1
    abs(iso - valp2) < eps && return p2
    abs(valp1-valp2) < eps && return p1
    mu = (iso - valp1) / (valp2 - valp1)
    p = p1 + mu * (p2 - p1)

    return p
end

"""
    mc_vert_points(xi,yi,zi, scale, origin, ::Type{VertType})

Returns a tuple of 8 points corresponding to each corner of a cube
"""
@inline function mc_vert_points(xi,yi,zi, scale, origin, ::Type{VertType}) where VertType
    (VertType(xi-1,yi-1,zi-1) .* scale .+ origin,
     VertType(xi  ,yi-1,zi-1) .* scale .+ origin,
     VertType(xi  ,yi  ,zi-1) .* scale .+ origin,
     VertType(xi-1,yi  ,zi-1) .* scale .+ origin,
     VertType(xi-1,yi-1,zi  ) .* scale .+ origin,
     VertType(xi  ,yi-1,zi  ) .* scale .+ origin,
     VertType(xi  ,yi  ,zi  ) .* scale .+ origin,
     VertType(xi-1,yi  ,zi  ) .* scale .+ origin)
end
