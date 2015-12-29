
# Linearly interpolate the position where an isosurface cuts
# an edge between two vertices, each with their own scalar value
function vertex_interp(iso,p1,p2,valp1,valp2, eps = 0.00001)

    abs(iso - valp1) < eps && return p1
    abs(iso - valp2) < eps && return p2
    abs(valp1-valp2) < eps && return p1
    mu = (iso - valp1) / (valp2 - valp1)
    p = p1 + mu * (p2 - p1)

    return p
end
