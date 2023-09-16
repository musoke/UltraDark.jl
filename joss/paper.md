---
title: 'UltraDark.jl: A Julia package for simulation of cosmological scalar fields'
tags:
  - Julia
  - cosmology
  - dark matter
  - inflation
  - scalar field
  - dynamics
authors:
  - name: Nathan Musoke
    orcid:  0000-0001-9839-9256
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: "Department of Physics and Astronomy, University of New Hampshire, USA"
   index: 1
date: 15 September 2023
bibliography: paper.bib

---

# Summary

`UltraDark.jl` is a Julia package for the simulation of cosmological scalar fields.
Scalar fields are proposed solutions to two of the fundamental questions in cosmology: the nature of dark matter and the universe's initial conditions.
Modeling their dynamics requires solving the Gross-Pitaevskii-Poisson equations, which is analytically challenging.
This makes simulations essential to understanding the dynamics of cosmological scalar fields.
`UltraDark.jl` is an open, performant and user friendly option for solving these equations numerically.


# Statement of need

Scalar fields are ubiquitous in physics, as solutions to partial differential equations describing the spatial variation of physical quantities.
As dark matter candidates, scalar fields including axion-like particles (ALPs) would explain the nature of the missing 85% of the universe's matter [@Planck:2015fie; @Matos:1999et; @Hu:2000ke; @Adams:2022pbo].
As inflaton candidates, scalar fields are proposed to cause a phase of accelerated expansion that sets the stage for big bang nucleosynthesis [@Guth:1980zm; @Linde:1981mu; @Albrecht:1982wi; @Amin:2011hj].
In each case, a scalar field $\psi(t, \mathbf{x})$ represents the density of particles as a function of space and time.


Subject to reasonable conditions, a cosmological scalar field $\psi$ whose particles have mass $m$ obeys the Gross-Pitaevskii equation for the scalar field,
\begin{equation}
	\label{eq:gpp}
	i \hbar \frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2 m a(t)^2} \nabla^2 \psi + m \Phi \psi
    \;,
\end{equation}
coupled to the Poisson equation for the gravitational potential $\Phi$
\begin{equation}
	\label{eq:poisson}
	\nabla^2 \Phi = \frac{4\pi G}{a(t)} m {|\psi|}^2
    \;,
\end{equation}
where $a(t)$ is the scale factor characterising the expansion of the universe.


These equations are difficult to solve analytically -- even static equilibrium solutions do not have a closed form -- and necessitate the use of computer simulations.
There are codes which solving \autoref{eq:gpp} and \autoref{eq:poisson} with different methods and in different domains, including `PyUltraLight` [@Edwards:2018ccc], a code written in Chapel [@padmanabhan2020simulating], `AxioNyx` [@Schwabe:2020eac], SCALAR [@Mina:2019ekb], i-SPin 2 [@Jain:2022agt; @Jain:2023qty]; see @Zhang:2018ghp for an overview.


`UltraDark.jl` solves \autoref{eq:gpp} and \autoref{eq:poisson} with a pseudo-spectral symmetrized split-step method, in which each time step consists of four sub-steps:
\begin{gather}
    \label{eq:update_phase_1}
    \psi \to \exp\left( - i\frac{h}{2} \Phi \right) \psi
    \\
    \label{eq:update_density}
    \psi \to \mathcal{F}^{-1}\left\{ \exp\left(-i h \frac{k^2}{2} \right) \mathcal{F}\left\{ \psi \right\} \right\}
    \\
    \Phi = \mathcal{F}^{-1}\left\{ - \frac{4\pi}{a k^2} \mathcal{F}\left\{ |\psi|^2\right\}\right\}
    \\
    \label{eq:update_phase_2}
    \psi \to \exp\left( - i\frac{h}{2} \Phi \right) \psi
    ,
\end{gather}
where $\mathcal{F}$ is a Fourier transform, $k$ are the corresponding frequencies, and $h$ is the time step.
`UltraDark.jl` has adaptive time steps which allow it to accelerate simulations while preserving numerical convergence.
This is particularly useful in an expanding universe, where the time step is roughly $h \propto a^2$.
Such time steps result in orders of magnitude speedups when simulating collapse of an inflaton field in the early universe [@Musoke:2019ima].


![Wall time for a single time step, as a function of number of CPUs. The points represent measured times and the lines represent theoretical $1/\text{\#CPU}$ scalings. The circles and solid line are for grids constructed from `Array`s and the triangles and dashed lines are for MPI-distributed `PencilArray`s.\label{fig:cpus}](../benchmarks/time_step/cpus.pdf){ width=100% }

Julia [@Julia-2017] has seen increasing use in scientific computing; see for example @Eschle:2023ikn and @JuliaBiologists for overviews of its use in high energy physics and biology.
The use of Julia is one of the choices that separates `UltraDark.jl` from similar codes.
`UltraDark.jl` uses Julia's has rich parallelism capabilities.
The `Threads.@threads` macro provides simple parallelisation of `for` loops.
Folds.jl enables simple parallelisation of reduction operations [@folds].
In a cluster environment, `PencilArrays.jl` and `PencilFFTs.jl` enable straightforward cross-node parallelism, a capability that is challenging to reproduce in Python [@PencilArrays; @PencilFFTs].
The scaling of this two approaches to parallelism is demonstrated in \autoref{fig:cpus}.


The features described above have allowed collaborators and I to produce results presented in 4 publications, each exploring the small scale structure of ultralight dark matter.
We have used `UltraDark.jl` to explore tidal disruption in dark matter halos comprised of self-interacting ALPs [@Glennon:2022huu], perform the first simulations of multi-species ALPs with intra- and inter-species interactions [@Glennon:2023jsp] and discover a novel mechanism for vortex stabilisation in scalar dark matter [@Glennon:2023oqa].
Work in preparation examines the effect of self-interactions on dynamical friction [@2023xxx].
More sample output can be found in @anim:vortex and @anim:multi.


# Acknowledgements

I thank Arka Banerjee, Ben Chandran, Richard Easther, Mateja Gosenca, Shaun Hotchkiss, Emily Kendall, Anthony Mirasola, Ethan Nadler, Mark Neyrinck, Chanda Prescod-Weinstein, Yourong Frank Wang, Risa Wechsler, Luna Zagorac for discussion and comments on various stages of the development of `UltraDark.jl`.
Particular thanks to Noah Glennon as the first user of `UltraDark.jl`.

Computations were performed on Marvin, a Cray CS500 supercomputer at UNH supported by the NSF MRI program under grant AGS-1919310.
The author wishes to acknowledge the use of New Zealand eScience Infrastructure (NeSI) high performance computing facilities, consulting support and/or training services as part of this research. New Zealand's national facilities are provided by NeSI and funded jointly by NeSI's collaborator institutions and through the Ministry of Business, Innovation & Employment's Research Infrastructure programme.
This work was performed in part at Aspen Center for Physics, which is supported by National Science Foundation under Grant No. PHY-1607611.
This work was partially supported by a grant from the Sloan Foundation.


# References
