---
title: 'Gyselalib++: A Portable C++ Library for Semi-Lagrangian Gyrokinetic Simulations'
tags:

- C++
- HPC
- plasma physics
- gyrokinetics

authors:

- given-names: Emily
  surname: Bourne
  affiliation: 1
  orcid: 0000-0002-3469-2338

- given-names: Virginie
  surname: Grandgirard
  affiliation: 2
  orcid: 0000-0001-7821-9107

- given-names: Yuuichi
  surname: Asahi
  affiliation: 3
  orcid: 0000-0002-9997-1274

- given-names: Julien
  surname: Bigot
  affiliation: 3
  orcid: 0000-0002-0015-4304

- given-names: Peter
  surname: Donnel
  affiliation: 2
  orcid: 0000-0002-6669-416X

- given-names: Alexander
  surname: Hoffmann
  affiliation: 4
  orcid: 0009-0005-0736-6122

- given-names: Abdelhadi
  surname: Kara
  affiliation: 2
  orcid: 0009-0000-9693-5828

- given-names: Philipp
  surname: Krah
  affiliation: 2
  orcid: 0000-0001-8982-4230

- given-names: Etienne
  surname: Malaboeuf
  affiliation: 5
  orcid: 0009-0007-3320-1406

- given-names: Kevin
  surname: Obrejan
  affiliation: 2
  orcid: 0000-0002-1906-4181

- given-names: Thomas
  surname: Padioleau
  affiliation: 3
  orcid: 0000-0001-5496-0013

- given-names: Matthieu
  surname: Protais
  affiliation: 2
  orcid: 0009-0003-5852-7446

- given-names: Pauline
  surname: Vidal
  affiliation: 4
  orcid: 0009-0008-7233-6351

affiliations:

- name: SCITAS, EPFL, CH-1015 Lausanne, Switzerland
  index: 1
- name: CEA, IRFM, 13108 Saint-Paul-lez-Durance Cedex, France
  index: 2
- name: >-
    Université Paris-Saclay, UVSQ, CNRS, CEA,
    Maison de la Simulation,
    91191, Gif-sur-Yvette, France
  index: 3
- name: >-
    Max-Planck-Institut für Plasmaphysik, Garching,
    Germany
  index: 4
- name: CINES, France
  index: 5

date: 24 April 2025
bibliography: paper.bib

---

## Summary

Gyselalib++ provides the mathematical building blocks to construct gyrokinetic plasma simulations in C++, simulating a distribution function discretised in phase space on a fixed grid.
It relies on the DDC library [@ddc] to statically type the discretisation dimensions; thus preventing many common sources of errors.
Via DDC, Gyselalib++ also leverages the Kokkos framework [@trott2022], ensuring portability across both CPU and various GPU architectures.
The library provides a variety of tools including semi-Lagrangian advection operators, solvers for partial differential equations (PDEs), quadrature rules, and time steppers.
All operators are designed to function correctly or raise compiler errors on non-orthonormal coordinate systems thanks to the static typing.

## Statement of Need

Plasma simulations are essential for the development of magnetic confinement fusion devices for energy production.
The low collisionality of such plasmas make kinetic models a judicious choice.
In particular gyrokinetic theory [@brizard2007; @krommes2012], which reduces the 6D problem to a 5D problem by removing high frequency dynamics, is a popular tool for such simulations [@garbet2010].
Despite the reduction in dimensionality such simulations still require massively powerful high-performance computing (HPC) resources.
For ITER-sized simulations, exascale resources would still be required.

The pre-existing Gysela code [@grandgirard2016], written in Fortran, originally aimed to simulate plasma in the core region of a tokamak using semi-Lagrangian advection with a distribution function discretised in phase space on a fixed grid.
This approach was shown to work well and allowed the study of many interesting physical phenomena [@Esteve2018;@Sarazin2021;@DifPradalier2022].
However expanding this code to use more complex mathematical methods such as non-uniform points (vital for handling the different magnitudes of physical quantities in the core and edge regions), and increasingly complex geometries (such as D-shape geometries, geometries including both closed and open field lines, and potentially stellarator geometries) has proved to be challenging and sometimes error-prone.
The challenges of such extensions are further amplified when trying to organise such a code for use on new GPU architectures, necessary for exascale simulations.
This is a challenge shared by other gyrokinetic codes [@trilaksono2023].

In the case of Gysela, several of the required changes would have affected a large percentage of the code base.
This means that the effort required was comparable to a complete rewrite, without the advantages such a rewrite can bring.
For example, we have been able to capitalise on C++'s strengths by using template programming to enforce the correctness of the implemented equations.
A common source of error is writing equations with implicit assumptions, such as assuming an orthonormal coordinate system, or specific properties like those of a circular coordinate system.
In Gyselalib++, equations are either expressed in tensor notation, so that they are either accurate for all geometries or do not compile, or they explicitly state their dependencies.
C++ further enables us to add static assertions for cases with restricted applicability to prevent their misuse.

In contrast to Gysela, Gyselalib++ has been conceived as a library, similar to the SeLaLib Fortran library [@selalib], whose independent elements are each unit tested and can be combined to build a final simulation.
This design makes the library more versatile, since its elements are not tied to a specific simulation and can be adapted to different needs.
The shared elements also provide more confidence in the reliability of the implementation as they can prove their validity across multiple applications.

Gyselalib++ includes a range of reusable mathematical operators for plasma simulations.
These include, but are not limited to, semi-Lagrangian advection schemes, numerical quadrature, differential operators (e.g. finite difference methods), and solvers for common partial differential equations.
A complete list of available operators can be found in the documentation[^doc].
Many of these tools are designed to work on a variety of grids, including non-uniform grids, which are especially important in edge-region simulations.
The library also supports MPI-based parallelism, either with distributed operators or with transpositions between different multi-rank storage layouts.

[^doc]: <https://gyselax.github.io/gyselalibxx/>

The VOICE code [@bourne2023] has already been rewritten in C++ using the mathematical tools provided by Gyselalib++.
Several common simulations including Landau damping (in 2D or 4D Cartesian phase-space coordinates), a bump-on-tail instability (in 2D Cartesian phase-space coordinates), and a guiding-centre model (on polar coordinates) have also been implemented.
While these examples are included primarily for illustration, they also serve as valuable test-beds for developing and validating new numerical methods.
As such, these examples are ideal for mathematicians looking to validate new methods in realistic, publication-ready test cases.

## Acknowledgements

This project has received funding from the EoCoE-III project (Energy-oriented Centre of Excellence for Exascale HPC applications).
This work has been carried out within the framework of the EUROfusion Consortium, funded by the European Union via the Euratom Research and Training Program (Grant Agreement No. 101052200 — EUROfusion).
Views and opinions expressed are those of the author(s) only and do not necessarily reflect those of the European Union or the European Commission.
Neither the European Union nor the European Commission can be held responsible for them.
Emily Bourne’s salary was paid for by the EUROfusion Advanced Computing Hub (ACH).
This project was provided with computer and storage resources by GENCI CINES thanks to the grant 2024-A0160502224 on the supercomputer  on the supercomputer Adastra’s the GENOA partition.

## References
