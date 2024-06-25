using Documenter
using DocumenterInterLinks
using Literate
using UltraDark

links = InterLinks(
    "Julia" => (
        "https://docs.julialang.org/en/v1/",
        joinpath(@__DIR__, "src", "inventories", "Julia.toml"),
    ),
    "MPI" => "https://juliaparallel.org/MPI.jl/dev/",
    "PencilArrays" => "https://jipolanco.github.io/PencilArrays.jl/dev/",
    "PencilFFTs" => "https://jipolanco.github.io/PencilFFTs.jl/dev/",
)

fallbacks = ExternalFallbacks()

# Use Literate.jl to convert *.jl examples to markdown files included in docs
const LITERATE_OUTPUT = joinpath(@__DIR__, "src", "examples")
const LITERATE_INPUT = joinpath(@__DIR__, "..", "examples")

for (root, _, files) in walkdir(LITERATE_INPUT), file in files
    # ignore non julia files
    splitext(file)[2] == ".jl" || continue

    # full path to a literate script
    ipath = joinpath(root, file)

    # generated output path
    opath = splitdir(replace(ipath, LITERATE_INPUT => LITERATE_OUTPUT))[1]

    # generate the Literate output
    Literate.script(ipath, opath)
    Literate.markdown(ipath, opath)
    Literate.notebook(ipath, opath)
end

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
            "Output Configuration"=>"man/output.md",
            "Simulation Configuration"=>"man/config.md",
            "Summary Statistics"=>"man/summary.md",
        ],
        "Examples" =>
            Any["PencilGrids"=>"examples/soliton_velocity.md", "2D"=>"examples/2d.md"],
        "API" => "api.md",
    ],
    repo = Remotes.GitHub("musoke", "UltraDark.jl"),
    sitename = "UltraDark.jl",
    authors = "Nathan Musoke <nathan.musoke@gmail.com>",
    plugins = [links, fallbacks],
)

deploydocs(; repo = "https://github.com/musoke/UltraDark.jl", devbranch = "main")
