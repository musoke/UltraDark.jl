# UltraDark

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://musoke.github.io/UltraDark.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://musoke.github.io/UltraDark.jl/dev)
[![Build Status](https://github.com/musoke/UltraDark.jl/workflows/CI/badge.svg)](https://github.com/musoke/UltraDark.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/musoke/UltraDark.jl?svg=true)](https://ci.appveyor.com/project/musoke/UltraDark-jl)
[![Codecov](https://codecov.io/gh/musoke/UltraDark.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/musoke/UltraDark.jl)

Simulations of cosmological scalar fields inspired by [PyUltraLight](https://github.com/auckland-cosmo/PyUltraLight).

UltraDark uses [PencilFFTs](https://jipolanco.github.io/PencilFFTs.jl) to do Fourier transforms in an MPI environment.


## Installation

You will need [Julia](https://julialang.org/).

To install and runs tests, open the Julia REPL and enter `pkg` mode by pressing
`]`.
```julia
pkg> dev https://github.com/musoke/UltraDark.jl

pkg> test UltraDark
```

Run a IJulia/Jupyter notebook with
```julia
pkg> add IJulia

julia> using IJulia

julia> notebook()
```
