var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/#Docstrings","page":"API","title":"Docstrings","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [UltraDark, UltraDark.Config, UltraDark.Output]","category":"page"},{"location":"api/#UltraDark.AbstractGrids","page":"API","title":"UltraDark.AbstractGrids","text":"AbstractGrids\n\nAbstract type for grids containing simulation data\n\n\n\n\n\n","category":"type"},{"location":"api/#UltraDark.Grids","page":"API","title":"UltraDark.Grids","text":"Grids(length, resol::Int)\nGrids(length_tuple, resol_tuple::Tuple{Int, Int, Int})\n\nstruct containing grids used in a simulation\n\nExamples\n\nCreate an empty grid with length length and resolution resol\n\njulia> using UltraDark\n\njulia> length = 1;\n\njulia> resol = 16;\n\njulia> Grids(length, resol);\n\njulia> size(ans.ψx)\n(16, 16, 16)\n\nCreate an empty length[1]xlength[2]xlength[3] grid with resolution resol[1]xresol[2]xresol[3].\n\njulia> using UltraDark\n\njulia> Grids((1.0, 1.0, 0.5), (64, 64, 32));\n\njulia> size(ans.ψx)\n(64, 64, 32)\n\n\n\n\n\n","category":"type"},{"location":"api/#UltraDark.PencilGrids","page":"API","title":"UltraDark.PencilGrids","text":"PencilGrids(length, resol)\nPencilGrids(length_tuple, resol_tuple::Tuple{Int, Int, Int})\n\nstruct containing grids used in a simulation\n\nEach grid is a PencilArray, allowing multiprocess FFTs. This comes with significant overhead so is only useful when running in a multi-node environment.\n\nExamples\n\nCreate an empty grid with length length and resolution resol.  Uses PencilFFTs to create PencilArrays.\n\njulia> using UltraDark\n\njulia> len = 1;\n\njulia> resol = 16;\n\njulia> PencilGrids(len, resol);\n\n\nCreate an empty length[1]xlength[2]xlength[3] grid with resolution resol[1]xresol[2]xresol[3].\n\njulia> using UltraDark\n\njulia> PencilGrids((1.0, 1.0, 0.5), (64, 64, 32));\n\n\n\n\n\n\n","category":"type"},{"location":"api/#UltraDark.E_grav-Tuple{Any, Any}","page":"API","title":"UltraDark.E_grav","text":"E_grav(grids)\nE_grav(grids, psi)\n\nGravitational potential energy\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.E_gravity_density-Tuple{Number, Real}","page":"API","title":"UltraDark.E_gravity_density","text":"E_gravity_density(psi, Phi)\n\nGravitational energy density of field psi in gravitational potential Phi\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.E_kq-Tuple{Any}","page":"API","title":"UltraDark.E_kq","text":"E_kq(grids)\nE_kq(grids, psi)\n\nSum of kinetic and quantum energies\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.E_total-Tuple{Any}","page":"API","title":"UltraDark.E_total","text":"E_total(grids; constants=nothing)\n\nTotal energy of the scalar field: the sum of the kinetic and quantum energies.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.actual_time_step-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.actual_time_step","text":"actual_time_step(max_time_step, time_interval, time_step_options)\n\nActual size and number of time steps that should be taken if the maximum is max_time_step. No more than time_step_options.update_period steps should be taken, and they should fit in time_interval.\n\nExamples\n\njulia> using UltraDark: actual_time_step, TimeStepOptions\n\njulia> actual_time_step(0.11, 1, TimeStepOptions())\n(0.1, 10)\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.add_external_potential!-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.add_external_potential!","text":"add_external_potential!(t, grids, constants)\n\nAdd a gravitational potential to the grid. By default this does nothing, but can be overridden in multiple dispatch.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.angular_momentum-Tuple{Any}","page":"API","title":"UltraDark.angular_momentum","text":"angular_momentum(grids)\nangular_momentum(grids, ψx, ρx)\n\nCalculate total angular momentum\n\nReturns\n\nL: AbstractArray with length 3\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.angular_momentum_density-Tuple{Any}","page":"API","title":"UltraDark.angular_momentum_density","text":"angular_momentum_density(grids)\nangular_momentum_density(grids, ψx, ρx)\n\nCalculate angular momentum density at each grid point\n\nReturns\n\nL: AbstractArray with dimensions 3 x resolx x resoly x resol_z\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.auxiliary_step!-NTuple{4, Any}","page":"API","title":"UltraDark.auxiliary_step!","text":"auxiliary_step!(Δt, grids, t, constants)\nauxiliary_step!(Δt, grids, t, constants, s; a = 1.0)\n\nDo an auxiliary inner step. By default this does nothing, but can be overridden in multiple dispatch.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.azimuthal_angle-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.azimuthal_angle","text":"polar_angle(grids, r0)\n\nCalculate the azimuthal angle in spherical or cylindrical coordinates\n\nThis is \\phi in conventional physics notation.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.dV-Tuple{Any}","page":"API","title":"UltraDark.dV","text":"dV(grids)\n\nCalculate the volume of each grid cell\n\nExamples\n\njulia> using UltraDark\n\njulia> box_length = 1.0;\n\njulia> resol = 16;\n\njulia> g = Grids(box_length, resol);\n\njulia> dV(g) * resol^3 == box_length^3\ntrue\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.evolve_to!-NTuple{5, Any}","page":"API","title":"UltraDark.evolve_to!","text":"Evolve grids forward from t_start to t_end\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.good_phase_diff-Tuple{Any, Any}","page":"API","title":"UltraDark.good_phase_diff","text":"bad_phase_diff(grids, config)\n\nCheck if the phase of the ψ field is in a trustable regime.\n\nReturns\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.inner_step!-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.inner_step!","text":"inner_step!(Δt, grids, constants; a=1.0)\ninner_step!(Δt, grids, constants, s; a=1.0)\n\nPerform the \"inner\" time step in the symmetrized split-step Fourier method.\n\nThis step applies the diffusion terms and updates the gravitational potential.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.k_vec-Tuple{Any, Any}","page":"API","title":"UltraDark.k_vec","text":"k_vec(lengths, resols)\n\nCalculate the Fourier frequencies of a box with side lengths lengths and resolutions resols\n\nExamples\n\njulia> using UltraDark: k_vec\n\njulia> kvec = k_vec((2π, 2π, 2π), (4, 4, 4));\n\njulia> kvec[1]\n4-element AbstractFFTs.Frequencies{Float64}:\n  0.0\n  1.0\n -2.0\n -1.0\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.mass-Tuple{Any, Any}","page":"API","title":"UltraDark.mass","text":"mass(grids)\nmass(grids, rho)\n\nCalculate total mass of a density field\n\nExamples\n\njulia> using UltraDark\n\njulia> g = Grids(1.0, 16);\n\njulia> g.ρx .= 0.0;\n\njulia> g.ρx[1, 1, 1] = 1.0;\n\njulia> UltraDark.mass(g) == 1.0 * (1.0 / 16)^3\ntrue\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.max_normed_phase_diff-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.max_normed_phase_diff","text":"normed_max_phase_grad(grids)\n\nCompute maximum phase gradient of a grid\n\nNormalised to ignore large gradients in regions with low density.  These tend to be anomalous.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.max_time_step-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.max_time_step","text":"max_time_step(grids, a, external_states)\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.max_time_step_external-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.max_time_step_external","text":"max_time_step_external(grids, a, state)\n\nCalculate the maximum time step implied by an external state\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.max_time_step_grids-Tuple{Any, Any}","page":"API","title":"UltraDark.max_time_step_grids","text":"max_time_step_grids(grids, a)\nmax_time_step_grids(grids::PencilGrids, a)\n\nCalculate an upper bound on the time step from grid properties\n\nThis time step depends on the gravitational potential and the resolution.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.outer_step!-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.outer_step!","text":"outer_step!(Δt, grids, constants; a=1.0)\nouter_step!(Δt, grids, constants, s; a=1.0)\n\nPerform the \"outer\" time step in the symmetrized split-step Fourier method.\n\nThis step only updates the phase of ψ applying accelerations due to gravity, the amplitude is not changed.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.phase_diff-Tuple{Any, Any}","page":"API","title":"UltraDark.phase_diff","text":"phase_diff(field, dir)\n\nCompute point-to-point difference of phase on a grid along a direction\n\nReturns an array of size (size(field)[1], size(field)[2], size(field)[3]) containing differences in direction dir.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.phase_diff-Tuple{PencilArrays.PencilArray, Any}","page":"API","title":"UltraDark.phase_diff","text":"phase_diff(field::PencilArray, dir)\n\nCompute point-to-point difference of phase on a grid along a direction\n\nReturns an array of size (size(field)[1], size(field)[2], size(field)[3]) containing differences in direction dir.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.polar_angle-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.polar_angle","text":"polar_angle(grids, r0)\n\nCalculate the polar angle in spherical coordinates\n\nThis is \\theta in conventional physics notation.\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.radius_cylindrical-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.radius_cylindrical","text":"radius_cylindrical(grids, r0)\n\nCalculate the radial coordinate in cylindrical coordinates\n\nExamples\n\njulia> using UltraDark\n\njulia> import UltraDark: radius_cylindrical, azimuthal_angle\n\njulia> box_length = 1.0;\n\njulia> resol = 16;\n\njulia> g = Grids(box_length, resol);\n\njulia> all(radius_cylindrical(g) .* cos.(azimuthal_angle(g)) .≈ g.x)\ntrue\n\njulia> all(\n           radius_cylindrical(g, (0.0, 0.0, 1.0)) .* cos.(azimuthal_angle(g, (0.0, 0.0, 1.0))) .≈\n           g.x,\n       )\ntrue\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.radius_spherical-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.radius_spherical","text":"radius_spherical(grids, r0)\n\nCalculate the radial coordinate in a spherical coordinate system\n\nExamples\n\njulia> using UltraDark\n\njulia> import UltraDark: radius_spherical, polar_angle, azimuthal_angle\n\njulia> box_length = 1.0;\n\njulia> resol = 16;\n\njulia> g = Grids(box_length, resol);\n\njulia> all(radius_spherical(g) .* sin.(polar_angle(g)) .* cos.(azimuthal_angle(g)) .≈ g.x)\ntrue\n\njulia> all(\n           radius_spherical(g, (1.0, 0.0, 0.0)) .* sin.(polar_angle(g, (1.0, 0.0, 0.0))) .*\n           cos.(azimuthal_angle(g, (1.0, 0.0, 0.0))) .+ 1.0 .≈ g.x,\n       )\ntrue\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.take_steps!-NTuple{8, Any}","page":"API","title":"UltraDark.take_steps!","text":"Take n steps with time step Δt\n\nExamples\n\njulia> using UltraDark: take_steps!, Grids, OutputConfig, Config\n\njulia> take_steps!(\n           Grids(1.0, 16),\n           0,\n           0.5,\n           10,\n           OutputConfig(mktempdir(), []),\n           Config.constant_scale_factor,\n           nothing,\n           (),\n       )\n5.0\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.update_gravitational_potential!-Tuple{Any}","page":"API","title":"UltraDark.update_gravitational_potential!","text":"update_gravitational_potential!(grids; a = 1.0)\nupdate_gravitational_potential!(grids, constants; a = 1.0)\n\nUpdate density grids.ρx and the gravitational potential grids.Φx based on grids.ψx\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.Config.TimeStepOptions","page":"API","title":"UltraDark.Config.TimeStepOptions","text":"TimeStepOptions(update_period=10, multiplier=1.0)\n\nstruct containing options controlling the size and calculation of time steps.\n\nupdate_period::Int64 controls how many steps are taken before the timestep is updated.\n\nmultiplier::Float64 multiplies the calculated maximum time step by a constant\n\nSee also: [SimulationConfig](@\n\nExamples\n\njulia> using UltraDark\n\njulia> TimeStepOptions()\nTimeStepOptions(10, 1.0)\n\n\n\n\n\n","category":"type"},{"location":"api/#UltraDark.Output.OutputConfig","page":"API","title":"UltraDark.Output.OutputConfig","text":"OutputConfig\n\nstruct containing information about what to output.\n\nsummary_statistics should be another struct whose constructor generates summary statistics from t, a, Δt and a grids object. If it is to be used with a PencilGrids object, each field must be concrete and have a binary representation that MPI can handle.\n\n\n\n\n\n","category":"type"},{"location":"api/#UltraDark.Output.output_external_state-NTuple{4, Any}","page":"API","title":"UltraDark.Output.output_external_state","text":"output_external_state(external_state, output_config, step, index)\n\nOutput states other than the ψ field.\n\nBy default this does nothing, but can be overloaded.\n\nArguments\n\nindex::Integer index of state in external states\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.Output.output_external_state_header-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.Output.output_external_state_header","text":"output_external_state_header(state, output_config)\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.Output.output_external_states_headers-Tuple{Any, Any}","page":"API","title":"UltraDark.Output.output_external_states_headers","text":"output_external_states_headers(external_states, output_config)\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.Output.output_grids-Tuple{Any, Any, Any}","page":"API","title":"UltraDark.Output.output_grids","text":"output_grids(grids, output_config, step)\n\nWrite output from grids as specified in output_config\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.Output.output_output_times-Tuple{Any, Any}","page":"API","title":"UltraDark.Output.output_output_times","text":"output_output_times(output_times, output_config)\n\nOutput the times corresponding to slices\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.Output.output_state-NTuple{4, Any}","page":"API","title":"UltraDark.Output.output_state","text":"output_state(grids, external_states, output_config, step)\n\nWrite out the grids and possible external states\n\n\n\n\n\n","category":"method"},{"location":"api/#UltraDark.Output.output_xyz-Tuple{Any, Any}","page":"API","title":"UltraDark.Output.output_xyz","text":"output_xyz(grids, output_config)\n\nOutput the spatial coordinates defining the grid\n\n\n\n\n\n","category":"method"},{"location":"man/init/#Soliton-initialisation","page":"Soliton Initialisation","title":"Soliton initialisation","text":"","category":"section"},{"location":"man/init/","page":"Soliton Initialisation","title":"Soliton Initialisation","text":"Modules = [UltraDark.Initialise]","category":"page"},{"location":"man/init/#UltraDark.Initialise.add_fdm_soliton!-Tuple{UltraDark.AbstractGrids, Vararg{Any, 5}}","page":"Soliton Initialisation","title":"UltraDark.Initialise.add_fdm_soliton!","text":"add_fdm_soliton!(grids::AbstractGrids, [psi], mass, position, velocity, phase, t0)\n\nAdd a fuzzy dark matter soliton to grids.  The density profile is rescaled to the desired mass and the phase is set to the desired velocity.\n\nIf the argument psi is passed, the soliton is added to the array-like psi rather than than grids.ψx. grids is still used to calculate coordinates.\n\nNote that due to coarse grid effects, the mass of the added soliton may not match the input mass exactly.\n\nThe included density profile from comes from PyUltraLight.\n\nExamples\n\nThe following  script adds a soliton to a grid and checks that the resulting total mass is correct.\n\nusing UltraDark\nusing UltraDark.Initialise\n\nresol = 64\ngrids = Grids(10.0, resol)\n\nmass = 10.0\nposition = [0.0, 0.0, 0.0]\nvelocity = [0.0, 0.0, 0.0]\nphase = 0.0\nt0 = 0.0\n\nadd_fdm_soliton!(grids, mass, position, velocity, phase, t0)\n\n# Ensure density is correct\nUltraDark.update_gravitational_potential!(grids, ())\n\nactual_mass = UltraDark.mass(grids)\n\nisapprox(mass, actual_mass, rtol = 1e-3)\n\n# output\n\ntrue\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#Summary-statistics","page":"Summary Statistics","title":"Summary statistics","text":"","category":"section"},{"location":"man/summary/","page":"Summary Statistics","title":"Summary Statistics","text":"It is often useful to collect summary statistics as a simulation runs and avoid storing complete simulation boxes. The UltraDark.Summary module makes this easy.","category":"page"},{"location":"man/summary/#Example","page":"Summary Statistics","title":"Example","text":"","category":"section"},{"location":"man/summary/","page":"Summary Statistics","title":"Summary Statistics","text":"Say you want to extract the minimum density at each time step. Define a struct to contain this data:","category":"page"},{"location":"man/summary/","page":"Summary Statistics","title":"Summary Statistics","text":"struct MinDensity\n    ρx_min::Float64\nend","category":"page"},{"location":"man/summary/","page":"Summary Statistics","title":"Summary Statistics","text":"and a function to construct it","category":"page"},{"location":"man/summary/","page":"Summary Statistics","title":"Summary Statistics","text":"function MaxDensity(sim_time, a, Δt, grids, constants, external_states)\n    ρx_min = minimum(grids.ρx)\n\n    MinDensity(ρx_max)\nend","category":"page"},{"location":"man/summary/","page":"Summary Statistics","title":"Summary Statistics","text":"Passing this object into a simulation will add a column to $output_path/summary.csv containing the single value contained by MinDensity.","category":"page"},{"location":"man/summary/","page":"Summary Statistics","title":"Summary Statistics","text":"More complex statistics are also possible.","category":"page"},{"location":"man/summary/#Docstrings","page":"Summary Statistics","title":"Docstrings","text":"","category":"section"},{"location":"man/summary/","page":"Summary Statistics","title":"Summary Statistics","text":"Modules = [UltraDark.Summary]","category":"page"},{"location":"man/summary/#UltraDark.Summary","page":"Summary Statistics","title":"UltraDark.Summary","text":"Summary\n\nModule containing utilities for computing and outputing summary statistics at each time step.\n\n\n\n\n\n","category":"module"},{"location":"man/summary/#UltraDark.Summary.EnergyGravity","page":"Summary Statistics","title":"UltraDark.Summary.EnergyGravity","text":"EnergyGravity\n\nGravitational potential energy\n\n\n\n\n\n","category":"type"},{"location":"man/summary/#UltraDark.Summary.EnergyKineticQuantum","page":"Summary Statistics","title":"UltraDark.Summary.EnergyKineticQuantum","text":"EnergyKineticQuantum\n\nGravitational potential energy\n\n\n\n\n\n","category":"type"},{"location":"man/summary/#UltraDark.Summary.MaxDensity","page":"Summary Statistics","title":"UltraDark.Summary.MaxDensity","text":"MaxDensity\n\n\n\n\n\n","category":"type"},{"location":"man/summary/#UltraDark.Summary.MaxDensityIndex","page":"Summary Statistics","title":"UltraDark.Summary.MaxDensityIndex","text":"MaxDensityIndex\n\nThis struct contains 4 useful pieces of information: the maxiumum value and three indices.  Each is output to a summary.\n\n\n\n\n\n","category":"type"},{"location":"man/summary/#UltraDark.Summary.MeanDensity","page":"Summary Statistics","title":"UltraDark.Summary.MeanDensity","text":"MeanDensity\n\nSummary statistic containing the mean density and the number of cells over which it was calculated.\n\nExamples\n\njulia> using UltraDark\n\njulia> g = Grids(1.0, 16);\n\njulia> Summary.MeanDensity(g)\nUltraDark.Summary.MeanDensity(0.0, 4096)\n\n\n\n\n\n","category":"type"},{"location":"man/summary/#UltraDark.Summary.RmsDensityContrast","page":"Summary Statistics","title":"UltraDark.Summary.RmsDensityContrast","text":"RmsDensityContrast\n\n\n\n\n\n","category":"type"},{"location":"man/summary/#UltraDark.Summary.TotalMass","page":"Summary Statistics","title":"UltraDark.Summary.TotalMass","text":"TotalMass\n\nTotal mass on a grid\n\nExamples\n\njulia> using UltraDark\n\njulia> g = Grids(1.0, 16);\n\njulia> Summary.TotalMass(0.0, 1.0, 1e-1, g, nothing, ())\nUltraDark.Summary.TotalMass(0.0)\n\n\n\n\n\n","category":"type"},{"location":"man/summary/#UltraDark.Summary.WallTime","page":"Summary Statistics","title":"UltraDark.Summary.WallTime","text":"WallTime\n\nThe current time in the real world.\n\nExamples\n\nWallTime created after another contains a later time.\n\njulia> using UltraDark\n\njulia> t1 = Summary.WallTime();\n\njulia> t2 = Summary.WallTime(0.0, 1.0, 1e-1, Grids(1.0, 16), nothing, ());\n\njulia> t1.date <= t2.date\ntrue\n\n\n\n\n\n","category":"type"},{"location":"man/summary/#UltraDark.Summary.column_title-Union{Tuple{Val{T}}, Tuple{T}} where T","page":"Summary Statistics","title":"UltraDark.Summary.column_title","text":"column_title(summary_struct)\n\nGet column title\n\nAssumes that the desired column title is the name of the first field.  This can be refined for other structs by defining a specialization for the given datatype.\n\nNote that the input argument is of type Val{T}, so that one can dispatch on T.\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.column_titles-Tuple{Any}","page":"Summary Statistics","title":"UltraDark.Summary.column_titles","text":"column_titles(summaries)\n\nGenerate column titles from an iterable of summaries\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.generate_summary_row-Tuple{Any}","page":"Summary Statistics","title":"UltraDark.Summary.generate_summary_row","text":"generate_summary_row(summaries)\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.get_relevant_data-Tuple{Any}","page":"Summary Statistics","title":"UltraDark.Summary.get_relevant_data","text":"get_relevant_data(summary_struct)\n\nFormat column entry as string\n\nBy default, assumes that the desired data is in the first field of the struct. This can be refined for other structs by defining a specialization for the given datatype and returning a comma separated string.\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.output_summary_header-Tuple{Any}","page":"Summary Statistics","title":"UltraDark.Summary.output_summary_header","text":"output_summary_header(output_config)\n\nWrite a header for a summary file\n\nThe header contains labels for each column of the summary CSV file. This function overwrites the current contents of the file.\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.output_summary_row-NTuple{7, Any}","page":"Summary Statistics","title":"UltraDark.Summary.output_summary_row","text":"output_summary_row(grids, output_config, t, a, Δt)\n\nWrite a new row to the summary file\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.pool_summarystat-Tuple{UltraDark.Summary.MaxDensity, UltraDark.Summary.MaxDensity}","page":"Summary Statistics","title":"UltraDark.Summary.pool_summarystat","text":"pool_summarystat(S1::MaxDensity, S2::MaxDensity)\n\nMPI reduction operator for max density\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.pool_summarystat-Tuple{UltraDark.Summary.MaxDensityIndex, UltraDark.Summary.MaxDensityIndex}","page":"Summary Statistics","title":"UltraDark.Summary.pool_summarystat","text":"pool_summarystat(S1::MaxDensity, S2::MaxDensity)\n\nMPI reduction operator for max density index\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.pool_summarystat-Tuple{UltraDark.Summary.MeanDensity, UltraDark.Summary.MeanDensity}","page":"Summary Statistics","title":"UltraDark.Summary.pool_summarystat","text":"pool_summarystat(S1::MeanDensity, S2::MeanDensity)\n\nMPI reduction operator for mean density\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.pool_summarystat-Tuple{UltraDark.Summary.RmsDensityContrast, UltraDark.Summary.RmsDensityContrast}","page":"Summary Statistics","title":"UltraDark.Summary.pool_summarystat","text":"pool_summarystat(S1::RmsDensityContrast, S2::RmsDensityContrast)\n\nMPI reduction operator for summary statistics.\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.pool_summarystat-Tuple{UltraDark.Summary.ScaleFactor, UltraDark.Summary.ScaleFactor}","page":"Summary Statistics","title":"UltraDark.Summary.pool_summarystat","text":"pool_summarystat(S1::ScaleFactor, S2::ScaleFactor)\n\nMPI reduction operator for scale factor\n\nCheck that scale factors are equal and return\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.pool_summarystat-Tuple{UltraDark.Summary.SimulationTime, UltraDark.Summary.SimulationTime}","page":"Summary Statistics","title":"UltraDark.Summary.pool_summarystat","text":"pool_summarystat(S1::SimulationTime, S2::SimulationTime)\n\nMPI reduction operator for simulation time\n\nCheck that times are equal and return\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.pool_summarystat-Tuple{UltraDark.Summary.TimeStep, UltraDark.Summary.TimeStep}","page":"Summary Statistics","title":"UltraDark.Summary.pool_summarystat","text":"pool_summarystat(S1::TimeStep, S2::TimeStep)\n\nMPI reduction operator for scale factor\n\nCheck that scale factors are equal and return\n\n\n\n\n\n","category":"method"},{"location":"man/summary/#UltraDark.Summary.pool_summarystat-Tuple{UltraDark.Summary.WallTime, UltraDark.Summary.WallTime}","page":"Summary Statistics","title":"UltraDark.Summary.pool_summarystat","text":"pool_summarystat(S1::WallTime, S2::WallTime)\n\nMPI reduction operator for wall time\n\nReturn the time of the first argument.\n\n\n\n\n\n","category":"method"},{"location":"man/install/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"man/install/","page":"Installation","title":"Installation","text":"First, install Julia.","category":"page"},{"location":"man/install/","page":"Installation","title":"Installation","text":"To install and run tests, open the Julia REPL and enter pkg mode by pressing ].","category":"page"},{"location":"man/install/","page":"Installation","title":"Installation","text":"$ julia\n\npkg> dev UltraDark.jl\n\npkg> test UltraDark","category":"page"},{"location":"man/install/","page":"Installation","title":"Installation","text":"This may take a while to run. There are a lot of packages to compile and some of the tests are of real simulations.","category":"page"},{"location":"man/install/","page":"Installation","title":"Installation","text":"If the tests succeed, you're all set!","category":"page"},{"location":"man/install/","page":"Installation","title":"Installation","text":"Next, check out the overview or to jump to running the example notebooks in IJulia.","category":"page"},{"location":"man/overview/#Overview","page":"Overview","title":"Overview","text":"","category":"section"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"This page is a brief overview of how to use UltraDark. Other examples go into more detail.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"UltraDark is built around UltraDark.AbstractGrids objects that contain coordinates and fields, and keep track of the relations between them. There are two types of grids included with UltraDark. UltraDark.Grids are built around regular Julia Core.Arrays. These are the grids used in most of the examples in this documentation. UltraDark.PencilGrids is built around PencilArrays.PencilArray and PencilFFTs. PencilGrids are useful for taking advantage of MPI parallelism when running in a cluster environment.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"warning: Warning\nNote that there is significant overhead involved in using PencilGrids. It is best to stick to Grids unless you are running jobs accross multiple nodes.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"You can define your own subtype of UltraDark.AbstractGrids if you wish to take advantage of other forms of parallelism or change the dynamics of the fields.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"Let's start by creating a Grids object.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"using UltraDark\n\nconst resol = 64\nconst box_length = 10.0\n\ngrids = Grids(box_length, resol);","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"grids is initially empty, as we can see by checking its mass","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"UltraDark.update_gravitational_potential!(grids) # ensure the density is up to date\nUltraDark.mass(grids)","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"Initial conditions are set by modifying the field ψx of grids. Let's use UltraDark.Initialise.add_fdm_soliton! to add a soliton with nonzero velocity to the center of grid.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"const mass = 10\nconst position0 = [0, 0, 0]\nconst velocity = [1, 0, 0]\nconst phase = 0\nconst t0 = 0\n\nUltraDark.Initialise.add_fdm_soliton!(grids, mass, position0, velocity, phase, t0)","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"Now, as expected, the mass on grids is nonzero","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"UltraDark.update_gravitational_potential!(grids) # ensure the density is up to date\nUltraDark.mass(grids)","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"note: Note\nThis does not exactly equal mass because we have used a coarse grid.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"Let's also check the location of the soliton.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"initial_indices = argmax(grids.ρx) # hide\n@assert initial_indices == CartesianIndex(Int(resol//2),  Int(resol//2), Int(resol//2)) # hide\nargmax(grids.ρx)","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"The indices of the maximum are at half resol in each dimension; this is the center of the box.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"After creating a AbstractGrids on which a simulation will happen and adding some matter to it, one must specify how the simulation should be carried out. The most important details are when and where to write output, and when the simulation should end.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"const output_times = 0:0.5:2\nconst output_dir = joinpath(mktempdir(), \"output\")\noutput_config = OutputConfig(output_dir, output_times; box = true)","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"See UltraDark.Output.OutputConfig for more details of how to configure the output.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"Now we are ready to run a simulation. Running this line will likely take some time, especially if UltraDark.jl has not been precompiled by your Julia installation.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"@time simulate!(grids, output_config)","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"We initialised the soliton with a nonzero velocity. Let's check if the soliton moved during the simulation.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"@assert argmax(grids.ρx) != initial_indices # hide\nargmax(grids.ρx)","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"We can also check this by loading output files from output_dir and plotting them.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"using CairoMakie\nusing LaTeXStrings\nusing NPZ\n\nfig = Figure()\n\nfor (i, t_index) in enumerate([1, length(output_times)])\n    rho= npzread(joinpath(output_dir, \"rho_$(t_index).npy\"))\n    ax = Axis(fig[1, i], title=L\"$t = %$(output_times[t_index])$\", aspect=DataAspect())\n    heatmap!(grids.x[:, 1, 1], grids.y[1, :, 1], rho[:, :, Int(resol//2)])\nend\nfig\nsave(\"overview.png\", fig); nothing # hide","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"(Image: )","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"After understanding this overview you can browse the example notebooks to see more complex simulations and analysis. Please open an issue if you run into problems or have feature requests.","category":"page"},{"location":"#UltraDark.jl","page":"Home","title":"UltraDark.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"UltraDark.jl is a tool for simulation of cosmological scalar fields, inspired by PyUltraLight and designed to be simple to use and extend.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It solves the non-interacting Gross-Pitaevskii equation","category":"page"},{"location":"","page":"Home","title":"Home","text":"    i fracpartial psipartial t\n    =\n    - frac12m a^2 nabla^2 psi\n    + m psi Phi","category":"page"},{"location":"","page":"Home","title":"Home","text":"for a scalar field psi with mass m, coupled to Poisson's equation for its gravitational potential Phi","category":"page"},{"location":"","page":"Home","title":"Home","text":"    nabla^2 Phi = frac4 pi Ga rho = 4 pi psi^2\n    ","category":"page"},{"location":"","page":"Home","title":"Home","text":"where a(t) is the scale factor of the universe and G is Newton's constant.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Such a scalar field describes scalar dark matter (SDM) in models including ultralight dark matter (ULDM), fuzzy dark matter (FDM), axion-like particles (ALPs) and the like. It also describes an inflaton field in the reheating epoch of the early universe.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please open an issue if you run into problems or have feature requests.","category":"page"},{"location":"","page":"Home","title":"Home","text":"If UltraDark contributes to your research, please cite it.","category":"page"},{"location":"#Academic-articles-that-have-used-UltraDark.jl","page":"Home","title":"Academic articles that have used UltraDark.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"N. Glennon, E. O. Nadler, N. Musoke, A. Banerjee, C. Prescod-Weinstein and R. H. Wechsler, Tidal disruption of solitons in self-interacting ultralight axion dark matter, Phys. Rev. D 105, no.12, 123540 (2022) doi:10.1103/PhysRevD.105.123540, [arXiv:2205.10336 [astro-ph.CO]].","category":"page"},{"location":"","page":"Home","title":"Home","text":"N. Glennon, A. E. Mirasola, N. Musoke, M. C. Neyrinck and C. Prescod-Weinstein, Scalar dark matter vortex stabilization with black holes, JCAP 07, 004 (2023) doi:10.1088/1475-7516/2023/07/004, [arXiv:2301.13220 [astro-ph.CO]].\n<center>\n<iframe width=\"750\" height=\"300\" src=\"https://www.youtube.com/embed/DYeL5UHQjdE\" title=\"Figure 2: Stable vortex-soliton, with initial perturbations\" frameborder=\"0\" allow=\"accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share\" referrerpolicy=\"strict-origin-when-cross-origin\" allowfullscreen></iframe>\n</center>","category":"page"},{"location":"","page":"Home","title":"Home","text":"N.~Glennon, N.~Musoke and C.~Prescod-Weinstein, Simulations of multifield ultralight axionlike dark matter, Phys. Rev. D 107, no.6, 063520 (2023) doi:10.1103/PhysRevD.107.063520,  [arXiv:2302.04302 [astro-ph.CO]].\n<center>\n<iframe width=\"750\" height=\"350\" src=\"https://www.youtube.com/embed/IENq5imeIzE\" title=\"Collisions between unbound solitons with 1- and 2- species of dark matter and equal phase\" frameborder=\"0\" allow=\"accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share\" referrerpolicy=\"strict-origin-when-cross-origin\" allowfullscreen></iframe>\n</center>","category":"page"},{"location":"","page":"Home","title":"Home","text":"N. Glennon, N. Musoke, E. O. Nadler, C. Prescod-Weinstein, and R. H. Wechsler, Dynamical friction in self-interacting ultralight dark matter, Phys. Rev. D 109, no.6, 6, 063501 (2024) doi:10.1103/PhysRevD.109.063501, [arXiv:2312.07684 [astro-ph.CO]].\n<center>\n<iframe width=\"750\" height=\"350\" src=\"https://www.youtube.com/embed/27LzOpsgXxs\" title=\"Dynamical friction in self-interacting ultralight axion dark matter, Figure 6\" frameborder=\"0\" allow=\"accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share\" referrerpolicy=\"strict-origin-when-cross-origin\" allowfullscreen></iframe>\n</center>","category":"page"}]
}
