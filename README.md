# UltraDark

[![DOI](https://zenodo.org/badge/265130304.svg)](https://zenodo.org/badge/latestdoi/265130304)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://musoke.github.io/UltraDark.jl/stable)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://musoke.github.io/UltraDark.jl/dev)
[![Build Status](https://github.com/musoke/UltraDark.jl/workflows/CI/badge.svg)](https://github.com/musoke/UltraDark.jl/actions)
[![Codecov](https://codecov.io/gh/musoke/UltraDark.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/musoke/UltraDark.jl)

Simulations of cosmological scalar fields inspired by [PyUltraLight](https://github.com/auckland-cosmo/PyUltraLight).

UltraDark uses [PencilFFTs](https://jipolanco.github.io/PencilFFTs.jl/stable/) to do Fourier transforms in an MPI environment.


## Installation

You will need [Julia](https://julialang.org/).

To install and run tests, open the Julia REPL and enter `pkg` mode by pressing
`]`.
```julia
pkg> dev https://github.com/musoke/UltraDark.jl

pkg> test UltraDark
```

Run a Jupyter notebook with
```julia
pkg> add IJulia

julia> using IJulia

julia> notebook()
```

The [documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://musoke.github.io/UltraDark.jl/stable) has details of how to do more.
