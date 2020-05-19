using Documenter, JultraDark

makedocs(;
    modules=[JultraDark],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/musoke/JultraDark.jl/blob/{commit}{path}#L{line}",
    sitename="JultraDark.jl",
    authors="Nathan Musoke <n.musoke@auckland.ac.nz>",
    assets=String[],
)

deploydocs(;
    repo="github.com/musoke/JultraDark.jl",
)
