using Documenter, Example

push!(LOAD_PATH, "../src/")
using Gismo

# List of subsection pages
SUBSECTION_PAGES = [
    "gsCore/index.md",
    "gsNurbs/index.md"
]

makedocs(
    sitename = "Gismo.jl",
    modules  = [Gismo],
    pages = [
        "Home" => "index.md",
        "Modules" => SUBSECTION_PAGES
    ]
)
