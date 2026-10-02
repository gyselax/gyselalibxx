# Code Design and Layout

This page describes the programming style used in Gyselalib++ and how the repository is organised.

## Understanding Functional Programming in Gyselalib++

There are multiple coding paradigms possible when writing code in C++. Two of the most common are object-oriented programming and functional programming. You may be familiar with object-oriented programming however Gyselalib++ uses functional programming.

### Key Difference: Functional vs. Object-Oriented Programming (OOP)

**Object-Oriented Programming (OOP)** is built around the concept of objects that encapsulate data and behaviour. Developers used to OOP often think in terms of classes, inheritance, and polymorphism. While these are powerful concepts, they are not central to how Gyselalib++ is structured. This is notably because inheritance, and polymorphism are often incompatible with GPU programming.

**Functional Programming (FP)** aims to make code more readable by making it more obvious where data is modified. This is done by separating data storage from operators. This way, data can never be modified inside its class without passing through an operator. In our code the functional programming paradigm was chosen as it allows us to write the code in a way which more closely resembles the equations.

### Operators and Data Structures

When using functional programming classes fit into one of two groups:

- **Data Structures**:
    - Serve as containers for data, often passed into operators for processing.
    - Avoid embedding behaviour directly into data structures; instead, let operators handle transformations.

- **Operators**:
    - Should implement an `operator()` method, allowing them to be used as callable objects.
    - Their internal variables should remain constant to align with functional programming principles.
    - Break problems into small, composable functions. Avoid the temptation to design large, monolithic classes.
    - Can be thought of as a function. Classes are used instead so the functions called by this function can be decided at compile time.

By separating operators (behaviour) from data structures (data), you ensure clarity, reusability, and better adherence to the functional programming style.

For example consider the case of a semi-Lagrangian advection. The operator is similar to a function which modifies a data structure containing the distribution function. By using a class the method used to interpolate the function (e.g. spline interpolation or Lagrange interpolation) can be specified without duplicating code.

## Navigating the Gyselalib++ Codebase

The Gyselalib++ codebase is organised into several folders. The main folders of interest are:

- **src/** : This folder contains the source code for the library. It contains folders grouping the code by subject. The `src/` folder and each sub-folder contains a `README.md` detailing its contents. It is simple to navigate within this documentation : [src](../../src/README.md)
    - **geometry...** : The `src/` folder contains sub-folders whose names begin with `geometry`. These sub-folders contain code which is specific to a given geometry. Each geometry is defined by the dimensions on which its equations are defined. In these folders you will also find sub-folders grouping the code by subject. Each `geometry...` folder will contain a folder called `geometry` containing a file `geometry.hpp`. This file contains type aliases which are useful for this geometry. This includes the definition of the classes representing the dimensions. If a `geometry.hpp` is included in a file then this file can only be used for that specific geometry. Therefore files in a `geometry...` folder cannot be used for other geometries even if the two appear to be compatible at first glance (e.g. files from the `geometryXVx` folder cannot be used for the `geometryXYVxVy` simulations).

- **tests/** : This folder contains the unit tests for the library. The folder structure is similar to that found in the `src/` folder.

- **simulations/** : This folder contains examples of simulations. In order to use the Gyselalib++ library for research applications we usually create a new repository which contains this repository as a sub-module, however there are still some simulations hosted in this repository. This is notably the case for simulations which have known properties that can be used as end-to-end tests. For example the Landau damping case.
