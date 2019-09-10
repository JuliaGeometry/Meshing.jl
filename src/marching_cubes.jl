
#Look up Table
include("lut/mc.jl")

"""
`marching_cubes(sdf::SignedDistanceField, [iso = 0.0,] [MT = HomogenousMesh{Point{3,Float64},Face{3,Int}}])`

Construct a `HomogenousMesh` from a `SignedDistanceField` using the
marching cubes algorithm. This method is faster than Marching Tetrahedra
and generates fewer vertices and faces (about 1/4 as many).
However it may generate non-manifold meshes, while Marching
Tetrahedra guarentees a manifold mesh.
"""
function marching_cubes(sdf::SignedDistanceField{3,ST,FT},
                        iso=0.0,
                        MT::Type{M}=SimpleMesh{Point{3,Float64},Face{3,Int}},
                        eps=0.00001, reduceverts=true) where {ST,FT,M<:AbstractMesh}
    VertType, FaceType = _determine_types(M, FT)
    marching_cubes(sdf, VertType, FaceType, iso, M, eps, reduceverts)
end

function marching_cubes(sdf::SignedDistanceField{3,ST,FT}, ::Type{VertType}, ::Type{FaceType},
                               iso=0.0,
                               MT::Type{M}=SimpleMesh{Point{3,Float64},Face{3,Int}},
                               eps=0.00001, reduceverts=true, insidepositive=Val(false)) where {ST,FT,M<:AbstractMesh, VertType, FaceType}
    nx, ny, nz = size(sdf)
    h = HyperRectangle(sdf)
    w = widths(h)
    orig = origin(HyperRectangle(sdf))

    # we subtract one from the length along each axis because
    # an NxNxN SDF has N-1 cells on each axis
    s = VertType(w[1]/(nx-1), w[2]/(ny-1), w[3]/(nz-1))

    # arrays for vertices and faces
    vts = VertType[]
    fcs = FaceType[]
    mt = max(nx,ny,nz)
    reduceverts == true && sizehint!(vts, mt*mt*5)
    reduceverts == false && sizehint!(vts, mt*mt*6)
    sizehint!(fcs, mt*mt*2)
    vertlist = Vector{VertType}(undef, 12)
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
        cubeindex = insidepositive === Val(true) ? _get_cubeindex_pos(iso_vals, iso) : _get_cubeindex(iso_vals, iso)

        # Cube is entirely in/out of the surface
        (cubeindex == 0x00 || cubeindex == 0xff) && continue

        points = (VertType(xi-1,yi-1,zi-1) .* s .+ orig,
                  VertType(xi,yi-1,zi-1) .* s .+ orig,
                  VertType(xi,yi,zi-1) .* s .+ orig,
                  VertType(xi-1,yi,zi-1) .* s .+ orig,
                  VertType(xi-1,yi-1,zi) .* s .+ orig,
                  VertType(xi,yi-1,zi) .* s .+ orig,
                  VertType(xi,yi,zi) .* s .+ orig,
                  VertType(xi-1,yi,zi) .* s .+ orig)

        # Find the vertices where the surface intersects the cube
        find_vertices_interp!(vertlist, points, iso_vals, cubeindex, iso, eps)

        # Create the triangle
        reduceverts == true && _mc_unique_triangles!(vts, fcs, vertlist, cubeindex, FaceType)
        reduceverts == false && _mc_create_triangles!(vts, fcs, vertlist, cubeindex, FaceType)
    end
    MT(vts,fcs)
end


function marching_cubes(f::Function,
                        bounds::HyperRectangle,
                        ::Type{VertType}, ::Type{FaceType},
                        samples::NTuple{3,Int}=(256,256,256),
                        iso=0.0,
                        MT::Type{M}=SimpleMesh{Point{3,Float64},Face{3,Int}},
                        eps=0.00001, reduceverts=true, insidepositive=Val(false)) where {M<:AbstractMesh, VertType, FaceType}

    FT = eltype(VertType)

    nx, ny, nz = samples[1], samples[2], samples[3]
    w = VertType(widths(bounds))
    orig = VertType(origin(bounds))

    # we subtract one from the length along each axis because
    # an NxNxN SDF has N-1 cells on each axis
    s = VertType(w[1]/(nx-1), w[2]/(ny-1), w[3]/(nz-1))

    # arrays for vertices and faces
    vts = VertType[]
    fcs = FaceType[]
    mt = max(nx,ny,nz)
    reduceverts == true && sizehint!(vts, mt*mt*5)
    reduceverts == false && sizehint!(vts, mt*mt*6)
    sizehint!(fcs, mt*mt*2)
    vertlist = Vector{VertType}(undef, 12)
    iso_vals = Vector{FT}(undef,8)
    @inbounds for xi = 1:nx-1, yi = 1:ny-1, zi = 1:nz-1

        points = (VertType(xi-1,yi-1,zi-1) .* s .+ orig,
                  VertType(xi,yi-1,zi-1) .* s .+ orig,
                  VertType(xi,yi,zi-1) .* s .+ orig,
                  VertType(xi-1,yi,zi-1) .* s .+ orig,
                  VertType(xi-1,yi-1,zi) .* s .+ orig,
                  VertType(xi,yi-1,zi) .* s .+ orig,
                  VertType(xi,yi,zi) .* s .+ orig,
                  VertType(xi-1,yi,zi) .* s .+ orig)

        if zi == 1
            for i = 1:8
                iso_vals[i] = f(points[i])
            end
        else
            iso_vals[1] = iso_vals[5]
            iso_vals[2] = iso_vals[6]
            iso_vals[3] = iso_vals[7]
            iso_vals[4] = iso_vals[8]
            iso_vals[5] = f(points[5])
            iso_vals[6] = f(points[6])
            iso_vals[7] = f(points[7])
            iso_vals[8] = f(points[8])
        end

        #Determine the index into the edge table which
        #tells us which vertices are inside of the surface
        cubeindex = insidepositive === Val(true) ? _get_cubeindex_pos(iso_vals, iso) : _get_cubeindex(iso_vals, iso)

        # Cube is entirely in/out of the surface
        (cubeindex == 0x00 || cubeindex == 0xff) && continue

        # Find the vertices where the surface intersects the cube
        # TODO this can use the underlying function to find the zero.
        # The underlying space is non-linear so there will be error otherwise
        find_vertices_interp!(vertlist, points, iso_vals, cubeindex, iso, eps)

        # Create the triangle
        reduceverts == true && _mc_unique_triangles!(vts, fcs, vertlist, cubeindex, FaceType)
        reduceverts == false && _mc_create_triangles!(vts, fcs, vertlist, cubeindex, FaceType)

    end
    MT(vts,fcs)
end


"""
    _mc_create_triangles!(vts, fcs, vertlist, cubeindex, FaceType)

Create triangles by only adding vertices within the voxel.
Each face does not share a reference to a vertex with another face.
"""
@inline function _mc_create_triangles!(vts, fcs, vertlist, cubeindex, FaceType)
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
@inline function _mc_unique_triangles!(vts, fcs, vertlist, cubeindex, FaceType)
    fct = length(vts)

    vert_to_add = _mc_verts[cubeindex]
    # Each vertex list will have atleast 3 elements so we can
    # add them to the list immediately
    push!(vts, vertlist[vert_to_add[1]],
               vertlist[vert_to_add[2]],
               vertlist[vert_to_add[3]])

    for i = 4:count_ones(edge_table[cubeindex])
        elt = vert_to_add[i]
        push!(vts, vertlist[elt])
    end
    offsets = _mc_connectivity[cubeindex]

    # There is atleast one face so we can push it immediately
    push!(fcs, FaceType(fct+offsets[3], fct+offsets[2], fct+offsets[1]))

    iszero(offsets[4]) && return
    push!(fcs, FaceType(fct+offsets[6], fct+offsets[5], fct+offsets[4]))

    iszero(offsets[7]) && return
    push!(fcs, FaceType(fct+offsets[9], fct+offsets[8], fct+offsets[7]))

    iszero(offsets[10]) && return
    push!(fcs, FaceType(fct+offsets[12], fct+offsets[11], fct+offsets[10]))

    iszero(offsets[13]) && return
    push!(fcs, FaceType(fct+offsets[15], fct+offsets[14], fct+offsets[13]))

end

"""
    find_vertices_interp!(vertlist, points, iso_vals, cubeindex, iso, eps)

Find the vertices where the surface intersects the cube
"""
@inline function find_vertices_interp!(vertlist, points, iso_vals, cubeindex, iso, eps)
    if !iszero(edge_table[cubeindex] & 0x001)
        vertlist[1] =
            vertex_interp(iso,points[1],points[2],iso_vals[1],iso_vals[2], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x002)
        vertlist[2] =
            vertex_interp(iso,points[2],points[3],iso_vals[2],iso_vals[3], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x004)
        vertlist[3] =
            vertex_interp(iso,points[3],points[4],iso_vals[3],iso_vals[4], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x008)
        vertlist[4] =
            vertex_interp(iso,points[4],points[1],iso_vals[4],iso_vals[1], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x010)
        vertlist[5] =
            vertex_interp(iso,points[5],points[6],iso_vals[5],iso_vals[6], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x020)
        vertlist[6] =
            vertex_interp(iso,points[6],points[7],iso_vals[6],iso_vals[7], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x040)
        vertlist[7] =
            vertex_interp(iso,points[7],points[8],iso_vals[7],iso_vals[8], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x080)
        vertlist[8] =
            vertex_interp(iso,points[8],points[5],iso_vals[8],iso_vals[5], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x100)
        vertlist[9] =
            vertex_interp(iso,points[1],points[5],iso_vals[1],iso_vals[5], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x200)
        vertlist[10] =
            vertex_interp(iso,points[2],points[6],iso_vals[2],iso_vals[6], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x400)
        vertlist[11] =
            vertex_interp(iso,points[3],points[7],iso_vals[3],iso_vals[7], eps)
    end
    if !iszero(edge_table[cubeindex] & 0x800)
        vertlist[12] =
            vertex_interp(iso,points[4],points[8],iso_vals[4],iso_vals[8], eps)
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



function (::Type{MT})(df::SignedDistanceField{3,ST,FT}, method::MarchingCubes)::MT where {MT <: AbstractMesh, ST, FT}
    VertType, FaceType = _determine_types(MT, FT)
    marching_cubes(df, VertType, FaceType, method.iso, MT, method.eps, method.reduceverts, Val(method.insidepositive))
end

function (::Type{MT})(f::Function, h::HyperRectangle, size::NTuple{3,Int}, method::MarchingCubes)::MT where {MT <: AbstractMesh}
    VertType, FaceType = _determine_types(MT)
    marching_cubes(f, h, VertType, FaceType, size, method.iso, MT, method.eps, method.reduceverts, Val(method.insidepositive))
end

function (::Type{MT})(f::Function, h::HyperRectangle, method::MarchingCubes; size::NTuple{3,Int}=(128,128,128))::MT where {MT <: AbstractMesh}
    VertType, FaceType = _determine_types(MT)
    marching_cubes(f, h, VertType, FaceType, size, method.iso, MT, method.eps, method.reduceverts, Val(method.insidepositive))
end