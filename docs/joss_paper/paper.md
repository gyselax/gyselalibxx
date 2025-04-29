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

- given-names: Dorian
  surname: Midou
  affiliation: 2
  orcid: 0009-0009-9593-8176

- given-names: Mathieu
  surname: Peybernes
  affiliation: 1
  orcid: 0009-0000-4011-9047

- given-names: Matthieu
  surname: Protais
  affiliation: 2
  orcid: 0009-0003-5852-7446

- given-names: Kevin
  surname: Obrejan
  affiliation: 2
  orcid: 0000-0002-1906-4181

- given-names: Thomas
  surname: Padioleau
  affiliation: 3
  orcid: 0000-0001-5496-0013

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
It relies on the DDC library[@ddc] to statically type the discretisation dimensions; thus preventing many common sources of errors.
Via DDC, Gyselalib++ also leverages the Kokkos framework, ensuring portability across both CPU and various GPU architectures.
The library provides a variety of tools including semi-Lagrangian advection operators, solvers for partial differential equations (PDEs), quadrature rules, and time steppers.
All operators are designed to function correctly or raise compiler errors on non-orthonormal coordinate systems thanks to the static typing.

## Statement of Need

Plasma simulations are essential for the development of magnetic confinement fusion devices for energy production.
The low collisionality of such plasmas make kinetic models a judicious choice.
In particular gyrokinetic theory[@brizard2007; @krommes2012], which reduces the 6D problem to a 5D problem by removing high frequency dynamics, is a popular tool for such simlations [@garbet2010].
Despite the reduction in dimensionality such simulations still require massively powerful high-performance computing (HPC) resources.

The most popular numerical methods for such simulations are particle-in-cell (PIC) methods and Eulerian methods, however both methods have inherent disadvantages.
PIC methods are limited by numerical noise which, without complex noise reduction techniques, is only slowly attenuated $1 / \sqrt{N}$ as resource usage is increased.
Eulerian methods avoid marker sampling noise by discretising the distribution function on a fixed grid, however explicit time integration leads to a Courant-Friedrichs-Lewy (CFL) stability condition which can severely limit the maximum possible time step.
An alternative approach is to combine these techniques using semi-Lagrangian advection.

The pre-existing Gysela code[@grandgirard2016], written in Fortran, originally aimed to simulate plasma in the core region of a tokamak using semi-Lagrangian advection with a distribution function discretised in phase space on a fixed grid.
This approach was shown to work well and allowed the study of many interesting physical phenomena[@TODO1;@TODO2;@ETC].
However expanding this code to use more complex mathematical methods such as non-uniform points (vital for handling the different magnitudes of physical quantities in the core and edge regions), and increasingly complex geometries (such as D-shape geometries, geometries including both closed and open field lines, and potentially stellarator geometries) has proved to be challenging and sometimes error-prone.
The challenges of such extensions are further amplified when trying to organise such a code for use on new GPU architectures.
This is a challenge shared by other gyrokinetic codes[@trilaksono2023].

In the case of Gysela, several of the required changes would have affected a large percentage of the code base.
This means that the effort required was comparable to a complete rewrite, without the advantages such a rewrite can bring.
Gyselalib++ is one of the first gyrokinetic codes to be written from scratch in C++.
This has allowed us to capitalise on C++'s strengths by using template programming to enforce the correctness of the implemented equations.
A common source of error is writing equations with implicit assumptions, such as assuming an orthonormal coordinate system, or specific properties like those of a circular coordinate system.
In Gyselalib++, equations are either expressed in tensor notation, so that they are either accurate for all geometries or do not compile, or they explicitly state their dependencies.
C++ further enables us to add static assertions for cases with restricted applicability to prevent their misuse.

In contrast to Gysela, Gyselalib++ has been conceived as a library whose independent elements are each unit tested and can be combined to build a final simulation.
This design makes the library more versatile, since its elements are not tied to a specific simulation and can be adapted to different needs.
The shared elements also provide more confidence in the reliability of the implementation as they can prove their validity across multiple applications.

The VOICE code[@bourne2023] has already been rewritten in C++ using the mathematical tools provided by Gyselalib++.

## Acknowledgements

This project has received funding from the EoCoE-III project (Energy-oriented Centre of Excellence for Exascale HPC applications).
This work has been carried out within the framework of the EUROfusion Consortium, funded by the European Union via the Euratom Research and Training Program (Grant Agreement No. 101052200 — EUROfusion).
Views and opinions expressed are those of the author(s) only and do not necessarily reflect those of the European Union or the European Commission.
Neither the European Union nor the European Commission can be held responsible for them.
Emily Bourne’s salary was paid for by the EUROfusion Advanced Computing Hub (ACH).

## References
