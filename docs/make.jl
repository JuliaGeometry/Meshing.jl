using Documenter
using Meshing

makedocs(
    sitename = "Meshing",
    format = Documenter.HTML(),
    modules = [Meshing],
    pages = ["API" => "api.md",
             "Examples" => "examples.md",
             "Internals" => "internals.md"]
)

deploydocs(
    repo = "github.com/JuliaGeometry/Meshing.jl.git",
)
