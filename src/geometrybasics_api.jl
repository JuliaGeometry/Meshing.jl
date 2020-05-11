function GB.Mesh(df::GT.SignedDistanceField{3,ST,FT},
              method::AbstractMeshingAlgorithm;
              pointtype=GB.Point{3,Float32},
              facetype=GB.GLTriangleFace) where {ST, FT}

    h = df.bounds
    vts, fcs = isosurface(df.data, method, pointtype, facetype, origin=pointtype(GB.origin(h)), widths=pointtype(GB.widths(h)))
    return GB.Mesh(vts, fcs)
end

function GB.Mesh(f::Function,
                 h::GB.Rect,
                 samples::NTuple{3,T},
                 method::AbstractMeshingAlgorithm;
                 pointtype=GB.Point{3,Float32},
                 facetype=GB.GLTriangleFace) where {T <: Integer}
    vts, fcs = isosurface(f, method, pointtype, facetype; samples=samples,
                          origin=pointtype(GB.origin(h)), widths=pointtype(GB.widths(h)))

    return GB.Mesh(vts, fcs)
end

function GB.Mesh(f::Function, h::GB.Rect, method::AbstractMeshingAlgorithm;
              samples::NTuple{3,T}=_DEFAULT_SAMPLES,
              pointtype=GB.Point{3,Float32},
              facetype=GB.GLTriangleFace) where {T <: Integer}

    return GB.Mesh(f, h, samples, method; pointtype=pointtype, facetype=facetype)
end

function GB.Mesh(volume::AbstractArray{T, 3},
                 method::AbstractMeshingAlgorithm;
                 pointtype=GB.Point{3,Float32},
                 facetype=GB.GLTriangleFace,
                 kwargs...) where {T}

    vts, fcs = isosurface(volume, method, pointtype, facetype; kwargs...)
    return GB.Mesh(vts, fcs)
end
