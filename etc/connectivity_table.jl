include("../src/lut/mc.jl")


connectivity_table = []
for i in eachindex(tri_table)
    a = zeros(Int, 15)
    for j = 1:15
        elt = findfirst(x -> x == tri_table[i][j], [_mc_verts[i]...])
        if elt == nothing || _mc_verts[i][elt] == 0
            a[j] = 0
        else
            a[j] = elt
        end
    end
    push!(connectivity_table, tuple(a...))
end

for elt in connectivity_table
    println(string(elt)*",")
end

