
"""
    _determine_types_gb(meshtype, fieldtype=Float64, facelen=3)
Given a subtype of AbstractMesh, determine the
type of vertex/point and face to use for internal computations.
Preference is given to the types specified by the Mesh call,
and will default to the `FieldType` for `SignedDistanceField`,
and Point{3,Float64}/TriangleFace{Int} for direct function sampling.
"""
function _determine_types_gb(pointtype, facetype, fieldtype=Float64, facelen=3)
    # determine the point and face types
    # preference is given to the Mesh types
    # followed by SDF if unspecified
    VertType = if pointtype isa Nothing
        GB.Point{3, fieldtype}
    else
        pointtype
    end

    FaceType = if facetype isa Nothing
        GB.NgonFace{facelen, Int}
    else
        facetype
    end

    return VertType, FaceType
end

function mesh(df::GT.SignedDistanceField{3,ST,FT}, method::AbstractMeshingAlgorithm;
              pointtype=nothing, facetype=nothing, ) where {ST, FT}

    vertex_eltype = promote_type(FT, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types_gb(pointtype, facetype, vertex_eltype,
                                          GB.default_face_length(method))
    h = df.bounds
    vts, fcs = isosurface(df.data, method, VertType, FaceType, origin=VertType(GB.origin(h)), widths=VertType(GB.widths(h)))
    return GB.Mesh(vts, fcs)
end

function GB.mesh(f::Function, h::GB.Rect, samples::NTuple{3,T}, method::AbstractMeshingAlgorithm;
              pointtype=nothing, facetype=nothing) where {T <: Integer}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types_gb(pointtype, facetype, vertex_eltype,
                                          default_face_length(method))
    vts, fcs = isosurface(f, method, VertType, FaceType; samples=samples,
                          origin=VertType(GB.origin(h)), widths=VertType(GB.widths(h)))
    return GB.Mesh(vts, fcs)
end

function GB.mesh(f::Function, h::GB.Rect, method::AbstractMeshingAlgorithm;
              samples::NTuple{3,T}=_DEFAULT_SAMPLES, pointtype=nothing,
              facetype=nothing) where {T <: Integer}
    return mesh(f, h, samples, method; pointtype=pointtype, facetype=facetype)
end

function GB.mesh(volume::AbstractArray{T, 3}, method::AbstractMeshingAlgorithm; pointtype=nothing, facetype=nothing, kwargs...) where {T}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types_gb(pointtype, facetype, vertex_eltype,
                                          default_face_length(method))
    vts, fcs = isosurface(volume, method, VertType, FaceType; kwargs...)
    return GB.Mesh(vts, fcs)
end
