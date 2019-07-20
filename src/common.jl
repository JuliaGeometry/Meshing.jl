
"""
given `iso_vals` and iso, return an 8 bit value corresping
to each corner of a cube. 0 indicates in the isosurface
and 1 indicates outside the surface
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