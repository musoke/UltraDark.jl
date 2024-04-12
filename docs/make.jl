using Documenter
using DocumenterInterLinks
using UltraDark

links = InterLinks(
    "Julia" => (
        "https://docs.julialang.org/en/v1/",
        joinpath(@__DIR__, "src", "inventories", "Julia.toml"),
    ),
    "PencilArrays" => "https://jipolanco.github.io/PencilArrays.jl/dev/",
    "PencilFFTs" => "https://jipolanco.github.io/PencilFFTs.jl/dev/",
)

fallbacks = ExternalFallbacks()

makedocs(;
    modules = [UltraDark],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://musoke.github.io/UltraDark.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Installation"=>"man/install.md",
            "Overview"=>"man/overview.md",
            "Soliton Initialisation"=>"man/init.md",
            "Summary Statistics"=>"man/summary.md",
        ],
        "API" => "api.md",
    ],
    repo = Remotes.GitHub("musoke", "UltraDark.jl"),
    sitename = "UltraDark.jl",
    authors = "Nathan Musoke <nathan.musoke@gmail.com>",
    plugins = [links, fallbacks],
)

deploydocs(; repo = "github.com/musoke/UltraDark.jl", devbranch = "main")
