# Overview

Welcome to Gyselalib++! This page lists the documentation that new users and developers should read, in the order in which we recommend reading it.
Most of these pages are not specific to new users, so they are found in other sections of the documentation. You will probably come back to them regularly as you work with the code.

1. **Install Gyselalib++** : Follow the [installation instructions](../../toolchains/README.md) to clone the repository, set up your environment, compile the code and run the tests.

2. **Understand how the code is designed** : Read [Getting Started with Gyselalib++](getting_started.md) to learn why Gyselalib++ uses a functional programming style and how the repository is organised.

3. **Learn how data is represented** : Read [DDC in Gyselalib++](../core_concepts/DDC_in_gyselalibxx.md). Almost every line of Gyselalib++ uses the types described on this page (`Coord`, `Idx`, `IdxStep`, `IdxRange`, `Field`, ...) so it is essential reading before working on the code. You will find it in the [Core Concepts](../core_concepts/index.md) section when you need to refer back to it.

4. **Build a complete simulation** : Work through the [Landau damping tutorial](landau_damping_tutorial.md) to see how the building blocks of the library are combined to create a simulation.

5. **Explore the library** : Browse the [Building Blocks](../../src/README.md) to find out which methods are already implemented, and the [Simulations](../../simulations/README.md) to see more examples of how they are used.

6. **Prepare your first contribution** : Before writing code, read the [Contributing](../contributing/index.md) section, in particular the [coding standard](../contributing/CODING_STANDARD.md) and the guide to [creating a branch](../how_to/Using_git.md#branches).

If something goes wrong along the way, the [Troubleshooting](../troubleshooting/index.md) section lists solutions to common problems, organised by the error message you see.
