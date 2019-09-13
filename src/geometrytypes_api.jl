function (::Type{MT})(df::SignedDistanceField{3,ST,FT}, method::AbstractMeshingAlgorithm)::MT where {MT <: AbstractMesh, ST, FT}
    vertex_eltype = promote_type(FT, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT, vertex_eltype, default_face_length(method))
    h = df.bounds
    vts, fcs = isosurface(df.data, method, VertType, FaceType, origin=VertType(origin(h)), widths=VertType(widths(h)))
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function, h::HyperRectangle, samples::NTuple{3,T}, method::AbstractMeshingAlgorithm)::MT where {MT <: AbstractMesh, T <: Integer}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype, default_face_length(method))
    vts, fcs = isosurface(f, method, VertType, FaceType, samples=samples, origin=VertType(origin(h)), widths=VertType(widths(h)))
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function, h::HyperRectangle, method::AbstractMeshingAlgorithm; samples::NTuple{3,T}=_DEFAULT_SAMPLES)::MT where {MT <: AbstractMesh, T <: Integer}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype, default_face_length(method))
    vts, fcs = isosurface(f, method, VertType, FaceType, samples=samples, origin=VertType(origin(h)), widths=VertType(widths(h)))
    MT(vts, fcs)::MT
end

function (::Type{MT})(volume::AbstractArray{T, 3}, method::AbstractMeshingAlgorithm; vargs...) where {MT <: AbstractMesh, T}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype, default_face_length(method))
    vts, fcs = isosurface(volume, method, VertType, FaceType, vargs...)
    MT(vts, fcs)::MT
end
