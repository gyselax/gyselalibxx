# Developer's FAQ

## What is the difference between Debug and Release mode?

A Release execution is supposed to be the fastest execution possible based on the assumption that the code is bug-free.
A Debug execution is there to make it possible to catch bugs, potentially at the cost of some more execution time.

## Should I use abort, assert, or static\_assert to raise an error?

`static_assert` ensures that an error is raised during the compilation. If it is possible to use `static_assert` then this option should always be preferred as it prevents bad code from being generated and incurs no costs in the resulting code. In reality it is not always possible to use `static_assert`. In this case either `assert` or `abort` must be used depending on the condition being tested.

A `static_assert` statement is written as:
```cpp
   static_assert(condition, "This is an explanation of what went wrong. "+
                    "It will be displayed in the compilation output.");
```

`assert` ensures that an error is raised in Debug mode while `abort` ensures that an error is raised in Release mode.
In order to choose which of the two should be used it is therefore important to ask yourself the following questions:
1. Can the condition be true in a bug-free code?
2. Is the condition dependent on user input?

If the answer to either of these questions is yes then it may be better to use `abort`, otherwise `assert` should be preferred to avoid incurring an unnecessary cost in Release mode.

On CPU an assertion is written as:
```cpp
#include <cassert>

...

   // An explanation can be added here
   assert(condition);
```

On GPU an assertion is written as:
```cpp
   // An explanation can be added here
   KOKKOS_ASSERT(condition);
```

On CPU an abort statement is written as:
```cpp
if (condition) {
    throw std::runtime_error("This is an explanation of what went wrong. "+
                "Ideally it also explains how the user input can be changed to avoid this problem.");
}
```

On GPU an abort statement is written as:
```cpp
if (condition) {
    Kokkos::abort("This is an explanation of what went wrong. "+
                "Ideally it also explains how the user input can be changed to avoid this problem.");
}
```

Note also that the way of annotating such errors depends on how they are raised. This reflects the fact that `assert` errors should only be seen by developers while `abort` errors will also be seen by users.
Assertions indicate the line and file where they originated so it is easy for developers to locate any associated comments. This also gives them  starting point for how to debug the error.
In contrast `abort` errors are written into the line that raises the error. This means that the error is printed to the command line when the executable fails. The message should include enough information for the user to understand what is happening so they have an idea of how to change their input to fix this problem.
