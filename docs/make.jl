using Documenter
using Meshing

makedocs(
    sitename = "Meshing",
    format = Documenter.HTML(),
   #modules = [Meshing],
    pages = ["Index" => "index.md",
             "API" => "api.md",
             "Examples" => "examples.md"]
)

deploydocs(
    repo = "github.com/JuliaGeometry/Meshing.jl.git",
)
