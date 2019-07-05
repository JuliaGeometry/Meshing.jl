
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
                               eps=0.00001) where {ST,FT,M<:AbstractMesh}
    nx, ny, nz = size(sdf)
    h = HyperRectangle(sdf)
    w = widths(h)
    orig = origin(HyperRectangle(sdf))

    # we subtract one from the length along each axis because
    # an NxNxN SDF has N-1 cells on each axis
    s = Point{3,Float64}(w[1]/(nx-1), w[2]/(ny-1), w[3]/(nz-1))

    # arrays for vertices and faces
    vts = Point{3,Float64}[]
    fcs = Face{3,Int}[]
    mt = max(nx,ny,nz)
    sizehint!(vts, mt*mt*6)
    sizehint!(fcs, mt*mt*2)
    vertlist = Vector{Point{3,Float64}}(undef, 12)
    @inbounds for xi = 1:nx-1, yi = 1:ny-1, zi = 1:nz-1

        points = (Point{3,Float64}(xi-1,yi-1,zi-1) .* s .+ orig,
                  Point{3,Float64}(xi,yi-1,zi-1) .* s .+ orig,
                  Point{3,Float64}(xi,yi,zi-1) .* s .+ orig,
                  Point{3,Float64}(xi-1,yi,zi-1) .* s .+ orig,
                  Point{3,Float64}(xi-1,yi-1,zi) .* s .+ orig,
                  Point{3,Float64}(xi,yi-1,zi) .* s .+ orig,
                  Point{3,Float64}(xi,yi,zi) .* s .+ orig,
                  Point{3,Float64}(xi-1,yi,zi) .* s .+ orig)
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
        cubeindex = _mc_cubeindex(iso_vals, iso)

        # Cube is entirely in/out of the surface
        (cubeindex == 0x00 || cubeindex == 0xff) && continue

        # Find the vertices where the surface intersects the cube
        find_vertices_interp!(vertlist, points, iso_vals, cubeindex, iso, eps)

        # Create the triangle
        _mc_create_triangles!(vts, fcs, vertlist, cubeindex)
    end
    MT(vts,fcs)
end


function marching_cubes(f::Function,
                        bounds::HyperRectangle,
                        samples::NTuple{3,Int}=(256,256,256),
                        iso=0.0,
                        MT::Type{M}=SimpleMesh{Point{3,Float64},Face{3,Int}},
                        eps=0.00001) where {ST,FT,M<:AbstractMesh}
    nx, ny, nz = samples[1], samples[2], samples[3]
    w = widths(bounds)
    orig = origin(bounds)

    # we subtract one from the length along each axis because
    # an NxNxN SDF has N-1 cells on each axis
    s = Point{3,Float64}(w[1]/(nx-1), w[2]/(ny-1), w[3]/(nz-1))

    # arrays for vertices and faces
    vts = Point{3,Float64}[]
    fcs = Face{3,Int}[]
    mt = max(nx,ny,nz)
    sizehint!(vts, mt*mt*6)
    sizehint!(fcs, mt*mt*2)
    vertlist = Vector{Point{3,Float64}}(undef, 12)
    iso_vals = Vector{Float64}(undef,8)
    @inbounds for xi = 1:nx-1, yi = 1:ny-1, zi = 1:nz-1

     points = (Point{3,Float64}(xi-1,yi-1,zi-1) .* s .+ orig,
               Point{3,Float64}(xi,yi-1,zi-1) .* s .+ orig,
               Point{3,Float64}(xi,yi,zi-1) .* s .+ orig,
               Point{3,Float64}(xi-1,yi,zi-1) .* s .+ orig,
               Point{3,Float64}(xi-1,yi-1,zi) .* s .+ orig,
               Point{3,Float64}(xi,yi-1,zi) .* s .+ orig,
               Point{3,Float64}(xi,yi,zi) .* s .+ orig,
               Point{3,Float64}(xi-1,yi,zi) .* s .+ orig)

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
        cubeindex = _mc_cubeindex(iso_vals, iso)

        # Cube is entirely in/out of the surface
        (cubeindex == 0x00 || cubeindex == 0xff) && continue

        # Find the vertices where the surface intersects the cube
        # TODO this can use the underlying function to find the zero.
        # The underlying space is non-linear so there will be error otherwise
        find_vertices_interp!(vertlist, points, iso_vals, cubeindex, iso, eps)

        # Create the triangle
        _mc_create_triangles!(vts, fcs, vertlist, cubeindex)
    end
    MT(vts,fcs)
end

@inline function _mc_create_triangles!(vts, fcs, vertlist, cubeindex)
    fct = length(vts) + 3

    push!(vts, vertlist[tri_table[cubeindex][1]],
               vertlist[tri_table[cubeindex][2]],
               vertlist[tri_table[cubeindex][3]])
    push!(fcs, Face{3,Int}(fct, fct-1, fct-2))

    iszero(tri_table[cubeindex][4]) && return
    fct += 3
    push!(vts, vertlist[tri_table[cubeindex][4]],
               vertlist[tri_table[cubeindex][5]],
               vertlist[tri_table[cubeindex][6]])
    push!(fcs, Face{3,Int}(fct, fct-1, fct-2))

    iszero(tri_table[cubeindex][7]) && return
    fct += 3
    push!(vts, vertlist[tri_table[cubeindex][7]],
               vertlist[tri_table[cubeindex][8]],
               vertlist[tri_table[cubeindex][9]])
    push!(fcs, Face{3,Int}(fct, fct-1, fct-2))

    iszero(tri_table[cubeindex][10]) && return
    fct += 3
    push!(vts, vertlist[tri_table[cubeindex][10]],
               vertlist[tri_table[cubeindex][11]],
               vertlist[tri_table[cubeindex][12]])
    push!(fcs, Face{3,Int}(fct, fct-1, fct-2))

    iszero(tri_table[cubeindex][13]) && return
    fct += 3
    push!(vts, vertlist[tri_table[cubeindex][13]],
               vertlist[tri_table[cubeindex][14]],
               vertlist[tri_table[cubeindex][15]])
    push!(fcs, Face{3,Int}(fct, fct-1, fct-2))
end

@inline function _mc_cubeindex(iso_vals, iso)
    cubeindex = iso_vals[1] < iso ? 0x01 : 0x00
    iso_vals[2] < iso && (cubeindex |= 0x02)
    iso_vals[3] < iso && (cubeindex |= 0x04)
    iso_vals[4] < iso && (cubeindex |= 0x08)
    iso_vals[5] < iso && (cubeindex |= 0x10)
    iso_vals[6] < iso && (cubeindex |= 0x20)
    iso_vals[7] < iso && (cubeindex |= 0x40)
    iso_vals[8] < iso && (cubeindex |= 0x80)
    cubeindex
end

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

# Linearly interpolate the position where an isosurface cuts
# an edge between two vertices, each with their own scalar value
function vertex_interp(iso, p1, p2, valp1, valp2, eps = 0.00001)

    abs(iso - valp1) < eps && return p1
    abs(iso - valp2) < eps && return p2
    abs(valp1-valp2) < eps && return p1
    mu = (iso - valp1) / (valp2 - valp1)
    p = p1 + mu * (p2 - p1)

    return p
end

struct MarchingCubes{T} <: AbstractMeshingAlgorithm
     iso::T
     eps::T
end

MarchingCubes(iso::T1=0.0, eps::T2=1e-3) where {T1, T2} = MarchingCubes{promote_type(T1, T2)}(iso, eps)

function (::Type{MT})(df::SignedDistanceField, method::MarchingCubes)::MT where {MT <: AbstractMesh}
     marching_cubes(df, method.iso, MT, method.eps)
end

function (::Type{MT})(f::Function, h::HyperRectangle, size::NTuple{3,Int}, method::MarchingCubes)::MT where {MT <: AbstractMesh}
     marching_cubes(f, h, size, method.iso, MT, method.eps)
end

function (::Type{MT})(f::Function, h::HyperRectangle, method::MarchingCubes; size::NTuple{3,Int}=(128,128,128))::MT where {MT <: AbstractMesh}
     marching_cubes(f, h, size, method.iso, MT, method.eps)
end