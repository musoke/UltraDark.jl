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
 - name: Department of Physics & Astronomy, University of New Hampshire, USA
   index: 1
date: 1 June 2023
bibliography: paper.bib

---

# Summary


# Statement of need


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @ fidgit.

For a quick reference, the following citation commands can be used:

- @Glennon:2023oqa  ->  "Author et al. (2001)"
- [@Glennon:2023jsp] -> "(Author et al., 2001)"
- [Edwards:2018ccc; @Glennon:2022huu] -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](../benchmarks/time_step/cpus.pdf)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](../benchmarks/time_step/resol.pdf){ width=20% }

# Acknowledgements


# References
