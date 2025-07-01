# How to debug problems in Gyselalib++ code

This is a quick guide for approaching and fixing simple bugs in Gyselalib++.
Gyselalib++ runs on both CPU and GPU, debugging works best if you first confirm whether a bug appears on CPU (where bugs are easier to diagnose) before diving into GPU-specific issues.
Tools specific to CPU or GPU can be used to help debugging. The use of these tools is not described here but some tools are recommended. More detailed information about their usage can be found online.

## CPU Debug Build

Always Start with a CPU Debug build.

The Debug build turns off compiler optimisations that can hide bugs or make backtraces confusing.
It also activates assertions that other devs have added to ensure you are using their tools as intended.

CPU builds tend to give clearer error messages than GPU builds which makes debugging simpler. Further standard debugging tools can be used for CPU builds.

In the case of an assertion failure, the line where the error is located should be indicated. Inspecting the code here should give you an idea of the cause of the error. Often there will be a comment next to the assertion in the code. If the function is being misused a backtrace can help locate the calling functions to track back to the code that you have added.

### Creating a debug build

It is important to note that simply adding `-DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-O0 -g"` to the cmake command is usually **not** sufficient to create a debug build. If you are following the normal installation/compilation procedure as described in [Gyselalib++ Installation](../../toolchains/README.md) then the toolchain will override the build type. Please ensure that you are using a toolchain which has `debug` in the name. If this is not the case and you do not see a suitable toolchain then you should modify one of the existing toolchains by changing the line `set(CMAKE_BUILD_TYPE Release)` to `set(CMAKE_BUILD_TYPE Debug)`.

### Finding the crashing line

If you want to identify the line which crashes `gdb` is a common interactive tool which can be used for this purpose.

```bash
gdb --args ./gysela_exe config_file.yml
run
backtrace
```

Many tutorials exist for a more involved use of `gdb` which allows the use of break points and printing of variables and their contents.

### Segmentation fault or other memory issues

A segmentation fault is usually an indication of memory misuse. Most cases of memory misuse should raise an assertion in a debug build. For those that don't `valgrind` is a common tool for analysing this kind of problem.

```bash
valgrind ./gysela_exe config_file.yml
```

One thing to bear in mind is that memory misuse does not always lead to an immediate crash. Instead it can simply corrupt the memory leading to a crash later in the execution. The valgrind errors should therefore be examined in order.

## GPU Debug Build

If the code runs as expected on CPU but crashes on GPU then you can already disregard the sections of the code which always run on CPU. This can help narrow down the problem.

Before analysing the code ensure that there are no compilation warnings.
The GPU compiler is less verbose than the CPU compiler but the majority of the warnings that it manages to raise are critical so make sure all warnings are eliminated before trying to run the program.

The most common cause of GPU errors comes from trying to access an object saved on CPU from the GPU or from trying to access an object saved on GPU from the CPU.
This is important when handling a Kokkos object or anything containing a Kokkos object. In particular care should be taken when accessing a `Field` or using a **coordinate transformation operator**.
In both cases static assertions are available to let you check at compile time if you are allowed to use a particular object.

For example if you want to access a `DFieldX` from GPU you can use:

```cpp
static_assert(Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace, typename DFieldX::memory_space>::accessible);
```

Remember a `Field` is defined on the GPU by default and a `Field` defined on the CPU should have a type like `host_t<Field>`.

For a coordinate transformation operator it is harder to deduce the accessibility directly from the type as analytical transformations such as `CircularToCartesian` can be used from both the CPU and the GPU.
To check if a coordinate transformation operator whose type is denoted `Mapping` is accessible from GPU you can use:

```cpp
static_assert(is_accessible_v<Kokkos::DefaultExecutionSpace, Mapping>);
```

Similarly to check if it is accessible from CPU you can use:

```cpp
static_assert(is_accessible_v<Kokkos::DefaultHostExecutionSpace, Mapping>);
```

An object is being accessed from the GPU if you are inside a function tagged with `KOKKOS_FUNCTION` or `KOKKOS_INLINE_FUNCTION`, or if you are inside a lambda function tagged with `KOKKOS_LAMBDA` or `KOKKOS_CLASS_LAMBDA`.

### Finding the crashing line

`gdb` is a CPU only tool but on a Nvidia system the tool `cuda-gdb` is a drop-in replacement that works on GPU.

### Memory issues

On a Nvidia system the tool `compute-sanitizer` can be used similarly to valgrind to weed out memory issues.
