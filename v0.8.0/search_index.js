var documenterSearchIndex = {"docs":
[{"location":"#UltraDark.jl-1","page":"Home","title":"UltraDark.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [UltraDark]","category":"page"},{"location":"#UltraDark.Grids","page":"Home","title":"UltraDark.Grids","text":"struct containing grids used in a simulation\n\nExamples\n\njulia> using UltraDark\n\njulia> len = 1;\n\njulia> resol = 16;\n\njulia> Grids(len, resol);\n\n\n\n\n\n\n","category":"type"},{"location":"#UltraDark.Grids-Tuple{Any, Int64}","page":"Home","title":"UltraDark.Grids","text":"Grids(length, resol::Int)\n\nConstructor for Grids\n\nCreate an empty grid with length length and resolution resol\n\nExamples\n\njulia> using UltraDark\n\njulia> Grids(1.0, 64);\n\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.Grids-Tuple{Any, Tuple{Int64, Int64, Int64}}","page":"Home","title":"UltraDark.Grids","text":"Grids(length_tuple, resol_tuple::Tuple{Int, Int, Int})\n\nConstructor for Grids\n\nCreate an empty length[1]xlength[2]xlength[3] grid with resolution resol[1]xresol[2]xresol[3].\n\nExamples\n\njulia> using UltraDark\n\njulia> Grids((1.0, 1.0, 0.5), (64, 64, 32));\n\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.PencilGrids","page":"Home","title":"UltraDark.PencilGrids","text":"struct containing grids used in a simulation\n\nExamples\n\njulia> using UltraDark\n\njulia> len = 1;\n\njulia> resol = 16;\n\njulia> PencilGrids(len, resol);\n\n\n\n\n\n\n","category":"type"},{"location":"#UltraDark.PencilGrids-Tuple{Any, Int64}","page":"Home","title":"UltraDark.PencilGrids","text":"PencilGrids(length, resol)\n\nConstructor for PencilGrids\n\nCreate an empty grid with length length and resolution resol.  Uses PencilFFTs to create PencilArrays.\n\nExamples\n\njulia> using UltraDark\n\njulia> PencilGrids(1.0, 64);\n\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.PencilGrids-Tuple{Any, Tuple{Int64, Int64, Int64}}","page":"Home","title":"UltraDark.PencilGrids","text":"PencilGrids(length_tuple, resol_tuple::Tuple{Int, Int, Int})\n\nConstructor for PencilGrids\n\nCreate an empty length[1]xlength[2]xlength[3] grid with resolution resol[1]xresol[2]xresol[3]. Each grid is a PencilArray, allowing multiprocess FFTs.\n\nExamples\n\njulia> using UltraDark\n\njulia> PencilGrids((1.0, 1.0, 0.5), (64, 64, 32));\n\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.dV-Tuple{Any}","page":"Home","title":"UltraDark.dV","text":"dV(grids)\n\nCalculate the volume of each grid cell\n\nExamples\n\njulia> using UltraDark\n\njulia> box_length = 1.0;\n\njulia> resol = 16;\n\njulia> g = Grids(box_length, resol);\n\njulia> dV(g) * resol^3 == box_length^3\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.E_grav-Tuple{Any, Any}","page":"Home","title":"UltraDark.E_grav","text":"E_grav(grids)\nE_grav(grids, psi)\n\nGravitational potential energy\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.E_kq-Tuple{Any}","page":"Home","title":"UltraDark.E_kq","text":"E_kq(grids)\nE_kq(grids, psi)\n\nSum of kinetic and quantum energies\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.E_total-Tuple{Any}","page":"Home","title":"UltraDark.E_total","text":"E_total(grids; constants=nothing)\n\nTotal energy of the scalar field: the sum of the kinetic and quantum energies.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.actual_time_step-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.actual_time_step","text":"actual_time_step(max_time_step, time_interval, time_step_options)\n\nActual size and number of time steps that should be taken if the maximum is max_time_step. No more than time_step_options.update_period steps should be taken, and they should fit in time_interval.\n\nExamples\n\njulia> using UltraDark: actual_time_step, TimeStepOptions\n\njulia> actual_time_step(0.11, 1, TimeStepOptions())\n(0.1, 10)\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.add_external_potential!-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.add_external_potential!","text":"add_external_potential!(t, grids, constants)\n\nAdd a gravitational potential to the grid. By default this does nothing, but can be overridden in multiple dispatch.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.angular_momentum-Tuple{Any}","page":"Home","title":"UltraDark.angular_momentum","text":"angular_momentum(grids)\nangular_momentum(grids, ψx, ρx)\n\nCalculate total angular momentum\n\nReturns\n\nL: AbstractArray with length 3\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.angular_momentum_density-Tuple{Any}","page":"Home","title":"UltraDark.angular_momentum_density","text":"angular_momentum_density(grids)\nangular_momentum_density(grids, ψx, ρx)\n\nCalculate angular momentum density at each grid point\n\nReturns\n\nL: AbstractArray with dimensions 3 x resolx x resoly x resol_z\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.auxiliary_step!-NTuple{4, Any}","page":"Home","title":"UltraDark.auxiliary_step!","text":"auxiliary_step!(Δt, grids, t, constants)\nauxiliary_step!(Δt, grids, t, constants, s; a = 1.0)\n\nDo an auxiliary inner step. By default this does nothing, but can be overridden in multiple dispatch.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.azimuthal_angle-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.azimuthal_angle","text":"polar_angle(grids, r0)\n\nCalculate the azimuthal angle in spherical or cylindrical coordinates\n\nThis is \\phi in conventional physics notation.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.evolve_to!-NTuple{5, Any}","page":"Home","title":"UltraDark.evolve_to!","text":"Evolve grids forward from t_start to t_end\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.good_phase_diff-Tuple{Any, Any}","page":"Home","title":"UltraDark.good_phase_diff","text":"bad_phase_diff(grids, config)\n\nCheck if the phase of the ψ field is in a trustable regime.\n\nReturns\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.inner_step!-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.inner_step!","text":"inner_step!(Δt, grids, constants; a=1.0)\ninner_step!(Δt, grids, constants, s; a=1.0)\n\nPerform the \"inner\" time step in the symmetrized split-step Fourier method.\n\nThis step applies the diffusion terms and updates the gravitational potential.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.k_vec-Tuple{Any, Any}","page":"Home","title":"UltraDark.k_vec","text":"k_vec(lengths, resols)\n\nCalculate the Fourier frequencies of a box with side lengths lengths and resolutions resols\n\nExamples\n\njulia> using UltraDark: k_vec\n\njulia> kvec = k_vec((2π, 2π, 2π), (4, 4, 4));\n\njulia> kvec[1]\n4-element AbstractFFTs.Frequencies{Float64}:\n  0.0\n  1.0\n -2.0\n -1.0\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.mass-Tuple{Any, Any}","page":"Home","title":"UltraDark.mass","text":"mass(grids)\nmass(grids, rho)\n\nCalculate total mass of a density field\n\nExamples\n\njulia> using UltraDark\n\njulia> g = Grids(1.0, 16);\n\njulia> g.ρx .= 0.0;\n\njulia> g.ρx[1, 1, 1] = 1.0;\n\njulia> UltraDark.mass(g) == 1.0 * (1.0 / 16)^3\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.max_normed_phase_diff-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.max_normed_phase_diff","text":"normed_max_phase_grad(grids)\n\nCompute maximum phase gradient of a grid\n\nNormalised to ignore large gradients in regions with low density.  These tend to be anomalous.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.max_time_step-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.max_time_step","text":"max_time_step(grids, a, external_states)\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.max_time_step_external-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.max_time_step_external","text":"max_time_step_external(grids, a, state)\n\nCalculate the maximum time step implied by an external state\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.max_time_step_grids-Tuple{Any, Any}","page":"Home","title":"UltraDark.max_time_step_grids","text":"max_time_step_grids(grids, a)\nmax_time_step_grids(grids::PencilGrids, a)\n\nCalculate an upper bound on the time step from grid properties\n\nThis time step depends on the gravitational potential and the resolution.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.outer_step!-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.outer_step!","text":"outer_step!(Δt, grids, constants; a=1.0)\nouter_step!(Δt, grids, constants, s; a=1.0)\n\nPerform the \"outer\" time step in the symmetrized split-step Fourier method.\n\nThis step only updates the phase of ψ applying accelerations due to gravity, the amplitude is not changed.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.phase_diff-Tuple{Any, Any}","page":"Home","title":"UltraDark.phase_diff","text":"phase_diff(field, dir)\n\nCompute point-to-point difference of phase on a grid along a direction\n\nReturns an array of size (size(field)[1], size(field)[2], size(field)[3]) containing differences in direction dir.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.phase_diff-Tuple{PencilArrays.PencilArray, Any}","page":"Home","title":"UltraDark.phase_diff","text":"phase_diff(field::PencilArray, dir)\n\nCompute point-to-point difference of phase on a grid along a direction\n\nReturns an array of size (size(field)[1], size(field)[2], size(field)[3]) containing differences in direction dir.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.polar_angle-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.polar_angle","text":"polar_angle(grids, r0)\n\nCalculate the polar angle in spherical coordinates\n\nThis is \\theta in conventional physics notation.\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.radius_cylindrical-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.radius_cylindrical","text":"radius_cylindrical(grids, r0)\n\nCalculate the radial coordinate in cylindrical coordinates\n\nExamples\n\njulia> using UltraDark\n\njulia> import UltraDark: radius_cylindrical, azimuthal_angle\n\njulia> box_length = 1.0;\n\njulia> resol = 16;\n\njulia> g = Grids(box_length, resol);\n\njulia> all(radius_cylindrical(g) .* cos.(azimuthal_angle(g)) .≈ g.x)\ntrue\n\njulia> all(\n           radius_cylindrical(g, (0.0, 0.0, 1.0)) .* cos.(azimuthal_angle(g, (0.0, 0.0, 1.0))) .≈\n           g.x,\n       )\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.radius_spherical-Tuple{Any, Any, Any}","page":"Home","title":"UltraDark.radius_spherical","text":"radius_spherical(grids, r0)\n\nCalculate the radial coordinate in a spherical coordinate system\n\nExamples\n\njulia> using UltraDark\n\njulia> import UltraDark: radius_spherical, polar_angle, azimuthal_angle\n\njulia> box_length = 1.0;\n\njulia> resol = 16;\n\njulia> g = Grids(box_length, resol);\n\njulia> all(radius_spherical(g) .* sin.(polar_angle(g)) .* cos.(azimuthal_angle(g)) .≈ g.x)\ntrue\n\njulia> all(\n           radius_spherical(g, (1.0, 0.0, 0.0)) .* sin.(polar_angle(g, (1.0, 0.0, 0.0))) .*\n           cos.(azimuthal_angle(g, (1.0, 0.0, 0.0))) .+ 1.0 .≈ g.x,\n       )\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#UltraDark.take_steps!-NTuple{8, Any}","page":"Home","title":"UltraDark.take_steps!","text":"Take n steps with time step Δt\n\nExamples\n\njulia> using UltraDark: take_steps!, Grids, OutputConfig, Config\n\njulia> take_steps!(\n           Grids(1.0, 16),\n           0,\n           0.5,\n           10,\n           OutputConfig(mktempdir(), []),\n           Config.constant_scale_factor,\n           nothing,\n           (),\n       )\n5.0\n\n\n\n\n\n","category":"method"}]
}
