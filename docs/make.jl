using Documenter, JultraDark

makedocs(;
    modules=[JultraDark],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://musoke.github.io/JultraDark.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/musoke/JultraDark.jl/blob/{commit}{path}#L{line}",
    sitename="JultraDark.jl",
    authors="Nathan Musoke <n.musoke@auckland.ac.nz>",
)

deploydocs(;
    repo="github.com/musoke/JultraDark.jl",
)
