using Documenter, UltraDark

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
    repo = "https://github.com/musoke/UltraDark.jl/blob/{commit}{path}#L{line}",
    sitename = "UltraDark.jl",
    authors = "Nathan Musoke <nathan.musoke@gmail.com>",
)

deploydocs(; repo = "github.com/musoke/UltraDark.jl", devbranch = "main")
