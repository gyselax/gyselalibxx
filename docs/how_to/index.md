# How-to Guides

The pages in this section give instructions for carrying out specific tasks.

## Building and testing

- [Compile the code](../../toolchains/README.md#compilation) : Build Gyselalib++ using one of the pre-made toolchains.
- [Run the tests](../../toolchains/README.md#running-tests) : Check that your installation is working as expected.

## Debugging

- [Create a debug build](Debugging_workflow.md#creating-a-debug-build) : Make sure that assertions are activated and compiler optimisations are turned off.
- [Find the line where the code crashes](Debugging_workflow.md#finding-the-crashing-line) : Use `gdb` (or `cuda-gdb` on GPU) to obtain a backtrace.
- [Investigate a segmentation fault](Debugging_workflow.md#segmentation-fault-or-other-memory-issues) : Use `valgrind` (or `compute-sanitizer` on GPU) to find memory misuse.
- [Debug code which only fails on GPU](Debugging_workflow.md#gpu-debug-build) : Check that objects are accessible from the execution space where they are used.

## Performance

- [Profile your code with Kokkos Tools](profiling.md#general-profiling-with-kokkos-tools) : Use the tools provided by Kokkos to measure performance.
- [Time individual kernels](profiling.md#using-simplekerneltimer) : Use the `SimpleKernelTimer` to find the most expensive parts of your code.

## Code quality

- [Run the static analysis](Development_tips.md#cppcheck) : Use `cppcheck` to check for common coding errors before submitting your code.

## Git

- [Create a branch](Using_git.md#branches) : Name and organise your branches for development.
- [Handle submodules](Using_git.md#submodules) : Solve common problems with the git submodules used by Gyselalib++.
