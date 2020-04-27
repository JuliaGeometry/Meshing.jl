import GeometryBasics: mesh

function mesh(df::SignedDistanceField{3,ST,FT}, method::AbstractMeshingAlgorithm;
              pointtype=nothing, facetype=nothing, ) where {ST, FT}

    vertex_eltype = promote_type(FT, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(pointtype, facetype, vertex_eltype,
                                          default_face_length(method))
    h = df.bounds
    vts, fcs = isosurface(df.data, method, VertType, FaceType, origin=VertType(origin(h)), widths=VertType(widths(h)))
    return Mesh(vts, fcs)
end

function mesh(f::Function, h::Rect, samples::NTuple{3,T}, method::AbstractMeshingAlgorithm;
              pointtype=nothing, facetype=nothing) where {T <: Integer}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(pointtype, facetype, vertex_eltype,
                                          default_face_length(method))
    vts, fcs = isosurface(f, method, VertType, FaceType; samples=samples,
                          origin=VertType(origin(h)), widths=VertType(widths(h)))
    return Mesh(vts, fcs)
end

function mesh(f::Function, h::Rect, method::AbstractMeshingAlgorithm;
              samples::NTuple{3,T}=_DEFAULT_SAMPLES, pointtype=nothing,
              facetype=nothing) where {T <: Integer}
    return mesh(f, h, samples, method; pointtype=pointtype, facetype=facetype)
end

function mesh(volume::AbstractArray{T, 3}, method::AbstractMeshingAlgorithm; pointtype=nothing, facetype=nothing, kwargs...) where {T}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(pointtype, facetype, vertex_eltype,
                                          default_face_length(method))
    vts, fcs = isosurface(volume, method, VertType, FaceType; kwargs...)
    return Mesh(vts, fcs)
end
