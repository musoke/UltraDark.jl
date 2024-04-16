# UltraDark.jl

UltraDark.jl is a tool for simulation of cosmological scalar fields, inspired by [PyUltraLight](https://github.com/auckland-cosmo/PyUltraLight) and designed to be simple to use and extend.

It solves the non-interacting Gross-Pitaevskii equation
```math
    i \frac{\partial \psi}{\partial t}
    =
    - \frac{1}{2m a^2} \nabla^2 \psi
    + m \psi \Phi
```
for a scalar field ``\psi`` with mass ``m``, coupled to Poisson's equation for its gravitational potential ``\Phi``
```math
    \nabla^2 \Phi = \frac{4 \pi G}{a} \rho = 4 \pi |\psi|^2
    ,
```
where `a(t)` is the scale factor of the universe and ``G`` is Newton's constant.


Such a scalar field describes scalar dark matter (SDM) in models including ultralight dark matter (ULDM), fuzzy dark matter (FDM), axion-like particles (ALPs) and the like.
It also describes an inflaton field in the reheating epoch of the early universe.


Please [open an issue](https://github.com/musoke/UltraDark.jl/issues/new) if you run into problems or have feature requests.

If UltraDark contributes to your research, please cite it.


#### Academic articles that have used UltraDark.jl

- N. Glennon, E. O. Nadler, N. Musoke, A. Banerjee, C. Prescod-Weinstein and R. H. Wechsler,
  *Tidal disruption of solitons in self-interacting ultralight axion dark matter*,
  Phys. Rev. D 105, no.12, 123540 (2022)
  [doi:10.1103/PhysRevD.105.123540](https://doi.org/10.1103/PhysRevD.105.123540),
  [[arXiv:2205.10336 [astro-ph.CO]]](https://arxiv.org/abs/2205.10336).


- N. Glennon, A. E. Mirasola, N. Musoke, M. C. Neyrinck and C. Prescod-Weinstein,
  *Scalar dark matter vortex stabilization with black holes*,
  JCAP 07, 004 (2023)
  [doi:10.1088/1475-7516/2023/07/004](https://doi.org/10.1088/1475-7516/2023/07/004),
  [[arXiv:2301.13220 [astro-ph.CO]]](https://arxiv.org/abs/2301.13220).
  ```@raw html
  <center>
  <iframe width="750" height="300" src="https://www.youtube.com/embed/DYeL5UHQjdE" title="Figure 2: Stable vortex-soliton, with initial perturbations" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
  </center>
  ```


- N. Glennon, N. Musoke and C. Prescod-Weinstein,
  *Simulations of multifield ultralight axionlike dark matter*,
  Phys. Rev. D 107, no.6, 063520 (2023)
  [doi:10.1103/PhysRevD.107.063520](https://doi.org/10.1103/PhysRevD.107.063520), 
  [[arXiv:2302.04302 [astro-ph.CO]]](https://arxiv.org/abs/2302.04302).
  ```@raw html
  <center>
  <iframe width="750" height="350" src="https://www.youtube.com/embed/IENq5imeIzE" title="Collisions between unbound solitons with 1- and 2- species of dark matter and equal phase" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
  </center>
  ```


- N. Glennon, N. Musoke, E. O. Nadler, C. Prescod-Weinstein, and R. H. Wechsler,
  *Dynamical friction in self-interacting ultralight dark matter*,
  Phys. Rev. D 109, no.6, 6, 063501 (2024)
  [doi:10.1103/PhysRevD.109.063501](https://doi.org/10.1103/PhysRevD.109.063501),
  [[arXiv:2312.07684 [astro-ph.CO]]](https://arXiv.org/abs/2312.07684).
  ```@raw html
  <center>
  <iframe width="750" height="350" src="https://www.youtube.com/embed/27LzOpsgXxs" title="Dynamical friction in self-interacting ultralight axion dark matter, Figure 6" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
  </center>
  ```
