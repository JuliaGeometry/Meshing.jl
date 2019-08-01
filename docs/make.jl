using Documenter
using Meshing

makedocs(
    sitename = "Meshing",
    format = Documenter.HTML(),
    modules = [Meshing],
    pages = ["Internals" => "internals.md"]
)

deploydocs(
    repo = "github.com/JuliaGeometry/Meshing.jl.git",
)
