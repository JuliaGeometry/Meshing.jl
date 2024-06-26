using BenchmarkTools

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

a = rand(8)

@benchmark _get_cubeindex(a, 0.5)

@inline function _get_cubeindex2(iso_vals, iso)
    cubeindex = 0
    for i in 1:8
        cubeindex |= (iso_vals[i] < iso ? 1 : 0) << (i - 1)
    end
    cubeindex
end

@benchmark _get_cubeindex2(a, 0.5)
