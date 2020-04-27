abstract type AbstractDistanceField end
abstract type AbstractUnsignedDistanceField <: AbstractDistanceField end
abstract type AbstractSignedDistanceField <: AbstractDistanceField end

"""
A `SignedDistanceField` is a uniform sampling of an implicit function.
The `bounds` field corresponds to the sampling space intervals on each axis.
The `data` field represents the value at each point whose exact location
can be rationalized from `bounds`.
The type is parameterized by:
* `N` - The dimensionality of the sampling space.
* `SpaceT` - the type of the space where we will uniformly sample.
* `FieldT` - the type resulting from evaluation of the implicit function.
Note that decoupling the space and field types is useful since geometry can
be formulated with integers and distances can be measured with floating points.
"""
mutable struct SignedDistanceField{N,SpaceT,FieldT} <: AbstractSignedDistanceField
    bounds::Rect{N,SpaceT}
    data::Array{FieldT,N}
end

Base.size(s::SignedDistanceField) = size(s.data)
Base.size(s::SignedDistanceField, i...) = size(s.data, i...)
Base.getindex(s::SignedDistanceField, i) = getindex(s.data,i)
Base.getindex(s::SignedDistanceField, i, j) = getindex(s.data,i, j)
Base.getindex(s::SignedDistanceField, i, j, k) = getindex(s.data,i, j,k)
Base.getindex(s::SignedDistanceField, i...) = getindex(s.data,i...)
GeometryBasics.Rect(s::SignedDistanceField) = s.bounds

"""
Construct a `SignedDistanceField` by sampling a function over the `bounds`
at the specified `resolution` (default = 0.1). Note that the sampling grid
must be regular,
so a new Rect will be generated for the SignedDistanceField that
may have larger maximum bounds than the input Rect. The default
Field type is Float64, but this can be changed with the `fieldT` argument.
"""
function SignedDistanceField(f::Function,
                             bounds::Rect{3, T},
                             resolution=0.1,
                             fieldT=Float64) where T
    x_min, y_min, z_min = minimum(bounds)
    x_max, y_max, z_max = maximum(bounds)

    x_rng, y_rng, z_rng = maximum(bounds) - minimum(bounds)

    nx = ceil(Int, x_rng/resolution)
    ny = ceil(Int, y_rng/resolution)
    nz = ceil(Int, z_rng/resolution)

    vol = Array{fieldT}(undef, nx+1, ny+1, nz+1)

    nb_max = Vec(x_min + resolution*nx,
                 y_min + resolution*ny,
                 z_min + resolution*nz)

    for i = 0:nx, j = 0:ny, k = 0:nz
        x = x_min + resolution*i
        y = y_min + resolution*j
        z = z_min + resolution*k
        @inbounds vol[i+1,j+1,k+1] = f(Vec{3,fieldT}(x,y,z))
    end

    nb_min = minimum(bounds)
    SignedDistanceField{3,T,fieldT}(Rect{3, T}(nb_min, nb_max-nb_min), vol)
end

function SignedDistanceField(f::Function,
                             bounds::Rect{2,T},
                             resolution=0.1,
                             fieldT=Float64) where T
    x_min, y_min = minimum(bounds)
    x_max, y_max = maximum(bounds)

    x_rng, y_rng = maximum(bounds) - minimum(bounds)

    nx = ceil(Int, x_rng/resolution)
    ny = ceil(Int, y_rng/resolution)

    vol = Array{fieldT}(undef, nx+1, ny+1)

    nb_max = Vec(x_min + resolution*nx,
                 y_min + resolution*ny)

    for i = 0:nx, j = 0:ny
        x = x_min + resolution*i
        y = y_min + resolution*j
        @inbounds vol[i+1,j+1] = f(Vec{2,fieldT}(x,y))
    end

    nb_min = minimum(bounds)
    SignedDistanceField{2,T,fieldT}(Rect{2, T}(nb_min, nb_max-nb_min), vol)
end
