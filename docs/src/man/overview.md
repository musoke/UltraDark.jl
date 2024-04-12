# Overview

This page is a brief overview of how to use UltraDark.
Other examples go into more detail.


UltraDark is built around [`UltraDark.AbstractGrids`](@ref) objects that contain coordinates and fields, and keep track of the relations between them.
There are two types of grids included with UltraDark.
[`UltraDark.Grids`](@ref) are built around regular Julia [`Core.Array`](@extref)s.
These are the grids used in most of the examples in this documentation.
[`UltraDark.PencilGrids`](@ref) is built around [`PencilArrays.PencilArray`](@extref) and [`PencilFFTs`](@extref).
`PencilGrids` are useful for taking advantage of MPI parallelism when running in a cluster environment.

!!! warning
    Note that there is significant overhead involved in using `PencilGrids`.
    It is best to stick to `Grids` unless you are running jobs on multiple nodes.

You can define your own subtype of `AbstractGrids` if you wish to take advantage of other forms of parallelism or change the dynamics of the fields.


Lets start by creating a `Grids` object.
```@example 1; continued=true
using UltraDark

const resol = 64
const box_length = 10.0

grids = Grids(box_length, resol);
```
`grids` is initially empty, as we can see by checking its mass
```@example 1; continued=false
UltraDark.update_gravitational_potential!(grids) # ensure the density is up to date
UltraDark.mass(grids)
```


Initial conditions are set by modifying the field `ψx` of an `AbstractGrids` object.
Let's add a soliton with nonzero velocity to the center of `grid`.
```@example 1; continued=false
const mass = 10
const position0 = [0, 0, 0]
const velocity = [1, 0, 0]
const phase = 0
const t0 = 0

UltraDark.Initialise.add_fdm_soliton!(grids, mass, position0, velocity, phase, t0)
```
See [`UltraDark.Initialise.add_fdm_soliton!`](@ref) for more details of adding solitons and other examples for other ways of setting initial conditions.


Now, as expected, the mass on the grids is nonzero
```@example 1; continued=false
UltraDark.update_gravitational_potential!(grids) # ensure the density is up to date
UltraDark.mass(grids)
```

!!! note
    This does not exactly equal `mass` because we have used a coarse grid.

Lets also check the location of the soliton.
```@example 1; continued=false
initial_indices = argmax(grids.ρx) # hide
@assert initial_indices == CartesianIndex(Int(resol//2),  Int(resol//2), Int(resol//2)) # hide
argmax(grids.ρx)
```
The indices of the maximum are at half `resol` in each dimension; this is the center of the box.


After creating an `AbstractGrids` on which a simulation will happen and adding some matter to it, one must specify how the simulation should be carried out.
The most important details are when and where to write output, and when the simulation should end.
```@example 1; continued=true
const output_times = 0:0.5:2

const output_dir = joinpath(mktempdir(), "output")
output_config = OutputConfig(output_dir, output_times; box = true)
```
See [`UltraDark.Output.OutputConfig`](@ref) for more details of how to configure the output.


Now we are ready to run a simulation.
Running this line will likely take some time, especially if UltraDark.jl has not been precompiled by your Julia installation.
```@example 1
@time simulate!(grids, output_config)
```

We initialised the soliton with a nonzero velocity `velocity`.
Lets check if the soliton moved during the simulation.
```@example 1; continued=false
@assert argmax(grids.ρx) != initial_indices # hide
argmax(grids.ρx)
```

We can also check this by loading output files from `output_dir` and plotting them.
```@example 1; continued=false
using CairoMakie
using LaTeXStrings
using NPZ

fig = Figure()

for (i, t_index) in enumerate([1, length(output_times)])
    rho= npzread(joinpath(output_dir, "rho_$(t_index).npy"))
    ax = Axis(fig[1, i], title=L"$t = %$(output_times[t_index])$", aspect=DataAspect())
    heatmap!(grids.x[:, 1, 1], grids.y[1, :, 1], rho[:, :, Int(resol//2)])
end
fig
save("overview.png", fig); nothing # hide
```
![](overview.png)

After understanding this overview you can browse the [example notebooks](https://github.com/musoke/UltraDark.jl/tree/main/examples) to see more complex simulations and analysis.
Please [open an issue](https://github.com/musoke/UltraDark.jl/issues/new) if you run into problems or have feature requests.
