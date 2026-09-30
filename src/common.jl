
"""
    _get_cubeindex(iso_vals, iso)

given `iso_vals` and iso, return an 8 bit value corresponding
to each corner of a cube. In each bit position,
0 indicates in the isosurface and 1 indicates outside the surface,
where the sign convention indicates negative inside the surface
"""
@inline function _get_cubeindex(iso_vals, iso)
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

"""
    no_triangles(cubeindex)

Called after `_get_cubeindex`. Determines if a voxel index has triangles.
"""
@inline function _no_triangles(cubeindex::UInt8)
    cubeindex == 0x00 || cubeindex == 0xff
end

"""
    smooth_sdf(sdf; sigma=0.7)

`sdf` blurred by a Gaussian of `sigma` voxels, edges clamped. Suppresses the banding a
grid-aligned level set gives the extracted surface. `sigma < 0.3` returns a copy.
"""
function smooth_sdf(sdf::AbstractArray{T,3}; sigma::Real=0.7) where {T}
    sigma < 0.3 && return copy(sdf)
    r = ceil(Int, 3sigma)
    w = [exp(-(k / sigma)^2 / 2) for k in -r:r]
    w = float(T).(w ./ sum(w))
    a, b = similar(sdf, float(T)), similar(sdf, float(T))
    blur!(a, sdf, w, 1)
    blur!(b, a, w, 2)
    blur!(a, b, w, 3)
end

function blur!(dst, src, w, axis)
    r = length(w) ÷ 2
    n = size(src, axis)
    @inbounds for I in CartesianIndices(src)
        acc = zero(eltype(dst))
        for k in -r:r
            J = CartesianIndex(ntuple(a -> a == axis ? clamp(I[a] + k, 1, n) : I[a], Val(3)))
            acc += w[k+r+1] * src[J]
        end
        dst[I] = acc
    end
    dst
end
