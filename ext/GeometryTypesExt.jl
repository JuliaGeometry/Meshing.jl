module GeometryTypesExt

using Meshing: AbstractMeshingAlgorithm, _DEFAULT_SAMPLES, isosurface, default_face_length
using GeometryTypes

const GT = GeometryTypes

"""
_determine_types(meshtype, fieldtype=Float64, facelen=3)

Given a subtype of AbstractMesh, determine the
type of vertex/point and face to use for internal computations.

Preference is given to the types specified by the Mesh call,
and will default to the `FieldType` for `SignedDistanceField`,
and Point{3,Float64}/Face{3,Int} for direct function sampling.
"""
function _determine_types(meshtype, fieldtype=Float64, facelen=3)
    # determine the point and face types
    # preference is given to the Mesh types
    # followed by SDF if unspecified
    if GT.vertextype(meshtype) !== Any
        VertType = GT.vertextype(meshtype)
    else
        VertType = GT.Point{3, fieldtype}
    end
    if GT.facetype(meshtype) !== Any
        FaceType = GT.facetype(meshtype)
    else
        FaceType = GT.Face{facelen, Int}
    end
    VertType, FaceType
end


function (::Type{MT})(df::GT.SignedDistanceField{3,ST,FT},
                      method::AbstractMeshingAlgorithm) where {MT <: GT.AbstractMesh, ST, FT}
    vertex_eltype = promote_type(FT, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT, vertex_eltype, default_face_length(method))
    h = df.bounds
    vts, fcs = isosurface(df.data, method, VertType, FaceType, origin=VertType(GT.origin(h)), widths=VertType(GT.widths(h)))
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function,
                      h::GT.HyperRectangle,
                      samples::NTuple{3,T},
                      method::AbstractMeshingAlgorithm)where {MT <: GT.AbstractMesh, T <: Integer}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype, default_face_length(method))
    vts, fcs = isosurface(f, method, VertType, FaceType, samples=samples, origin=VertType(GT.origin(h)), widths=VertType(GT.widths(h)))
    MT(vts, fcs)::MT
end

function (::Type{MT})(f::Function,
                      h::GT.HyperRectangle,
                      method::AbstractMeshingAlgorithm;
                      samples::NTuple{3,T}=_DEFAULT_SAMPLES) where {MT <: GT.AbstractMesh, T <: Integer}
    MT(f, h, samples, method)
end

function (::Type{MT})(volume::AbstractArray{T, 3},
                      method::AbstractMeshingAlgorithm;
                      kwargs...) where {MT <: GT.AbstractMesh, T}
    vertex_eltype = promote_type(T, typeof(method.iso), typeof(method.eps))
    VertType, FaceType = _determine_types(MT,vertex_eltype, default_face_length(method))
    vts, fcs = isosurface(volume, method, VertType, FaceType; kwargs...)
    MT(vts, fcs)::MT
end

end