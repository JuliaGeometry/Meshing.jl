
#
# MarchingCubes
#

function (::Type{MT})(df::SignedDistanceField{3,ST,FT}, method::MarchingCubes)::MT where {MT <: AbstractMesh, ST, FT}
    VertType, FaceType = _determine_types(MT, FT)
    h = df.bounds
    vts, fcs = isosurface(df.data, method, VertType, FaceType, origin=VertType(origin(h)), widths=VertType(widths(h)))
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function, h::HyperRectangle, samples::NTuple{3,T}, method::MarchingCubes)::MT where {MT <: AbstractMesh, T <: Integer}
    VertType, FaceType = _determine_types(MT)
    vts, fcs = isosurface(f, method, VertType, FaceType, samples=samples, origin=VertType(origin(h)), widths=VertType(widths(h)))
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function, h::HyperRectangle, method::MarchingCubes; samples::NTuple{3,T}=(128,128,128))::MT where {MT <: AbstractMesh, T <: Integer}
    VertType, FaceType = _determine_types(MT)
    vts, fcs = isosurface(f, method, VertType, FaceType, samples=samples, origin=VertType(origin(h)), widths=VertType(widths(h)))
    MT(vts, fcs)::MT
end

function (::Type{MT})(volume::AbstractArray{T, 3}, method::MarchingCubes; vargs...) where {MT <: AbstractMesh, T}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype)
    vts, fcs = isosurface(volume, method, VertType, FaceType, vargs...)
    MT(vts, fcs)::MT
end


#
# MarchingTetrahedra
#

function (::Type{MT})(sdf::SignedDistanceField{3,ST,FT}, method::MarchingTetrahedra) where {ST, FT, MT <: AbstractMesh}
    vertex_eltype = promote_type(FT, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype)
    h = sdf.bounds
    vts, fcs = isosurface(sdf.data, method, VertType, FaceType, origin=VertType(origin(h)), widths=VertType(widths(h)))
    MT(vts, fcs)::MT
end

function (::Type{MT})(volume::AbstractArray{T, 3}, method::MarchingTetrahedra; vargs...) where {MT <: AbstractMesh, T}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype)
    vts, fcs = isosurface(volume, method, VertType, FaceType, vargs...)
    MT(vts, fcs)::MT
end


#
# NaiveSurfaceNets
#

function (::Type{MT})(sdf::SignedDistanceField{3,ST,FT}, method::NaiveSurfaceNets) where {MT <: AbstractMesh, ST,FT}

    VertType, FaceType = _determine_types(MT, FT, 4)
    bounds = sdf.bounds
    vts, fcs = isosurface(sdf.data, method, VertType, FaceType, origin=origin(bounds), widths=widths(bounds))
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function, bounds::HyperRectangle, samples::NTuple{3,Int}, method::NaiveSurfaceNets) where {MT <: AbstractMesh}
    VertType, FaceType = _determine_types(MT, Float64, 4)
    vts, fcs = isosurface(f, method, VertType, FaceType, samples=samples, origin=origin(bounds), widths=widths(bounds))
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function, bounds::HyperRectangle, method::NaiveSurfaceNets;samples::NTuple{3,Int}=(128,128,128)) where {MT <: AbstractMesh}
    MT(f,bounds,samples,method)
end

function (::Type{MT})(volume::AbstractArray{T, 3}, method::NaiveSurfaceNets; vargs...) where {MT <: AbstractMesh, T}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype, 4)
    vts, fcs = isosurface(volume, method, VertType, FaceType, vargs...)
    MT(vts, fcs)::MT
end
