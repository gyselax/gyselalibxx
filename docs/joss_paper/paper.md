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



##
