using Documenter, UltraDark

makedocs(;
    modules = [UltraDark],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://musoke.github.io/UltraDark.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
    repo = "https://github.com/musoke/UltraDark.jl/blob/{commit}{path}#L{line}",
    sitename = "UltraDark.jl",
    authors = "Nathan Musoke <n.musoke@auckland.ac.nz>",
)

deploydocs(; repo = "github.com/musoke/UltraDark.jl")
