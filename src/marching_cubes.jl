
#Look up Table
include("lut/mc.jl")

#=
Voxel corner and edge indexing conventions

        Z
        |

        5------8------8  
       /|            /|      
      5 |           7 |      
     /  9          /  12     
    6------6------7   |      
    |   |         |   |      
    |   1------4--|---4  -- Y
    10 /          11 /       
    | 1           | 3        
    |/            |/        
    2------2------3    

  /
 X
=#

function isosurface(sdf::AbstractArray{T,3}, method::MarchingCubes, X=-1:1, Y=-1:1, Z=-1:1) where {T}
    vts, fcs, _ = mc_isosurface(sdf, method, X, Y, Z, Val(false))
    vts, fcs
end

"""
    isosurface_normals(sdf, method::MarchingCubes, X=-1:1, Y=-1:1, Z=-1:1) -> (vertices, faces, normals)

Like [`isosurface`](@ref), with a unit normal per vertex from the gradient of `sdf`.
"""
isosurface_normals(sdf::AbstractArray{T,3}, method::MarchingCubes, X=-1:1, Y=-1:1, Z=-1:1) where {T} =
    mc_isosurface(sdf, method, X, Y, Z, Val(true))

function mc_isosurface(sdf::AbstractArray{T,3}, method::MarchingCubes, X, Y, Z, ::Val{normals}) where {T,normals}
    nx, ny, nz = size(sdf)

    # find widest type
    FT = promote_type(eltype(first(X)), eltype(first(Y)), eltype(first(Z)), eltype(T), typeof(method.iso))

    vts = NTuple{3,float(FT)}[]
    fcs = NTuple{3,Int}[]
    nms = normals ? NTuple{3,float(FT)}[] : nothing

    xp = LinRange(first(X), last(X), nx)
    yp = LinRange(first(Y), last(Y), ny)
    zp = LinRange(first(Z), last(Z), nz)

    # the vertex of each crossed grid edge, when vertices are shared
    ids = method.reduceverts ? Dict{Int,Int}() : nothing
    mc_voxels!(vts, nms, fcs, ids, sdf, xp, yp, zp, method.iso)
    vts, fcs, nms
end

function mc_voxels!(vts, nms, fcs, ids, sdf, xp, yp, zp, iso)
    nx, ny, nz = size(sdf)
    h = (step(xp), step(yp), step(zp))
    idx = Vector{Int}(undef, 12)

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
        cubeindex = _get_cubeindex(iso_vals, iso)

        # Cube is entirely in/out of the surface
        _no_triangles(cubeindex) && continue

        points = mc_vert_points(xi, yi, zi, xp, yp, zp)
        grads = mc_vert_grads(nms, sdf, CartesianIndex(xi, yi, zi), h)

        # process the voxel
        mc_voxel!(vts, nms, fcs, ids, idx, cubeindex, (xi, yi, zi), (nx, ny), points, grads, iso, iso_vals)
    end
end

# no shared vertices and no normals: the original path
mc_voxel!(vts, ::Nothing, fcs, ::Nothing, idx, cubeindex, voxel, dims, points, grads, iso, iso_vals) =
    process_mc_voxel!(vts, fcs, cubeindex, points, iso, iso_vals)

function mc_voxel!(vts, nms, fcs, ids, idx, cubeindex, voxel, dims, points, grads, iso, iso_vals)
    @inbounds begin
        vert_to_add = _mc_verts[cubeindex]
        for i = 1:12
            vt = vert_to_add[i]
            iszero(vt) && break
            idx[i] = mc_vertex!(vts, nms, ids, vt, voxel, dims, points, grads, iso, iso_vals)
        end

        offsets = _mc_connectivity[_mc_eq_mapping[cubeindex]]
        push!(fcs, (idx[3], idx[2], idx[1]))
        for i in (1, 4, 7, 10)
            iszero(offsets[i]) && return
            push!(fcs, (idx[offsets[i+2]], idx[offsets[i+1]], idx[offsets[i]]))
        end
    end
end

# the index of edge `vt`'s vertex: a new one, or the one its grid edge already has
mc_vertex!(vts, nms, ::Nothing, vt, voxel, dims, points, grads, iso, iso_vals) =
    push_mc_vertex!(vts, nms, vt, points, grads, iso, iso_vals)
mc_vertex!(vts, nms, ids::Dict, vt, voxel, dims, points, grads, iso, iso_vals) =
    get!(() -> push_mc_vertex!(vts, nms, vt, points, grads, iso, iso_vals), ids, mc_edge_id(voxel, vt, dims))

function push_mc_vertex!(vts, nms, vt, points, grads, iso, iso_vals)
    a, b = _mc_edge_list[vt]
    push!(vts, vertex_interp(iso, points[a], points[b], iso_vals[a], iso_vals[b]))
    push_mc_normal!(nms, vertex_interp(iso, grads[a], grads[b], iso_vals[a], iso_vals[b]))
    length(vts)
end
function push_mc_vertex!(vts, nms, vt, points, ::Nothing, iso, iso_vals)
    a, b = _mc_edge_list[vt]
    push!(vts, vertex_interp(iso, points[a], points[b], iso_vals[a], iso_vals[b]))
    length(vts)
end

push_mc_normal!(nms, g) = push!(nms, all(iszero, g) ? (zero(g[1]), zero(g[1]), one(g[1])) : g ./ sqrt(sum(abs2, g)))

# each cube edge as (offset of its lower grid point, axis), so voxels sharing an edge give it one id
const mc_edge_offsets = ((0, 0, 0, 1), (1, 0, 0, 2), (0, 1, 0, 1), (0, 0, 0, 2),
                         (0, 0, 1, 1), (1, 0, 1, 2), (0, 1, 1, 1), (0, 0, 1, 2),
                         (0, 0, 0, 3), (1, 0, 0, 3), (1, 1, 0, 3), (0, 1, 0, 3))

@inline function mc_edge_id((xi, yi, zi), edge, (nx, ny))
    dx, dy, dz, axis = mc_edge_offsets[edge]
    axis + 3 * ((xi + dx - 1) + nx * ((yi + dy - 1) + ny * (zi + dz - 1)))
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

# corner offsets, in `mc_vert_points` order
const mc_corners = (CartesianIndex(0, 0, 0), CartesianIndex(1, 0, 0), CartesianIndex(1, 1, 0), CartesianIndex(0, 1, 0),
                    CartesianIndex(0, 0, 1), CartesianIndex(1, 0, 1), CartesianIndex(1, 1, 1), CartesianIndex(0, 1, 1))

mc_vert_grads(::Nothing, sdf, I, h) = nothing
mc_vert_grads(nms, sdf, I, h) = map(c -> sdf_grad(sdf, I + c, h), mc_corners)

sdf_grad(sdf, I, h) = ntuple(a -> sdf_deriv(sdf, I, a, h[a]), Val(3))

# 4th-order central difference, 2nd-order next to the border, one-sided on it
@inline function sdf_deriv(sdf, I, a, h)
    d = CartesianIndex(ntuple(b -> Int(b == a), Val(3)))
    i, n = I[a], size(sdf, a)
    @inbounds if 2 < i < n - 1
        (sdf[I-2d] - 8sdf[I-d] + 8sdf[I+d] - sdf[I+2d]) / 12h
    elseif 1 < i < n
        (sdf[I+d] - sdf[I-d]) / 2h
    elseif i == 1
        (sdf[I+d] - sdf[I]) / h
    else
        (sdf[I] - sdf[I-d]) / h
    end
end
