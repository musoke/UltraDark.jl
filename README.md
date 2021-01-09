# JultraDark

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://musoke.github.io/JultraDark.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://musoke.github.io/JultraDark.jl/dev)
[![Build Status](https://github.com/musoke/JultraDark.jl/workflows/CI/badge.svg)](https://github.com/musoke/JultraDark.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/musoke/JultraDark.jl?svg=true)](https://ci.appveyor.com/project/musoke/JultraDark-jl)
[![Codecov](https://codecov.io/gh/musoke/JultraDark.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/musoke/JultraDark.jl)
[![Build Status](https://api.cirrus-ci.com/github/musoke/JultraDark.jl.svg)](https://cirrus-ci.com/github/musoke/JultraDark.jl)

Simulations of cosmological scalar fields inspired by [PyUltraLight](https://github.com/auckland-cosmo/PyUltraLight).

Eventually may use [PencilFFTs](https://jipolanco.github.io/PencilFFTs.jl) to do Fourier transforms in distributed memory.


## Installation

You will need [Julia](https://julialang.org/).

To install and runs tests, open the Julia REPL and enter `pkg` mode by pressing
`]`.
```julia
pkg> dev https://github.com/musoke/JultraDark.jl

pkg> test JultraDark
```

Run a IJulia/Jupyter notebook with
```julia
pkg> add IJulia

julia> using IJulia

julia> notebook()
```
