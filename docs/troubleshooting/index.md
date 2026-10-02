# Troubleshooting

This page lists common problems, organised by the symptom or error message that you see.
If your problem is not listed here, please check the [issues](https://github.com/gyselax/gyselalibxx/issues) on GitHub or [report a bug](../CONTRIBUTING.md#bug-reports).

## Compilation errors

Search this table for a part of the error message reported by your compiler.
Compiler errors refer to types by their DDC names (e.g. `ddc::ChunkSpan` instead of `Field`). The [DDC and Gyselalib++ names](../core_concepts/ddc_names.md) page can be used to translate them.

| The error message contains | Solution |
| --- | --- |
| `The closure type for a lambda ... cannot be used in the template argument type of a __global__ function` | [Lambdas must capture by copy](Common_compilation_problems.md#the-closure-type-for-a-lambda-cannot-be-used-in-the-template-argument-type-of-a-__global__-function) |
| `Implicit capture of 'this' in extended lambda expression` (a warning which may lead to wrong results) | [Avoid capturing class members](Common_compilation_problems.md#implicit-capture-of-this-in-extended-lambda-expression) |
| `expression must be a modifiable lvalue` or `cannot be referenced -- it is a deleted function` | [Do not copy objects which own memory](Common_compilation_problems.md#accessing-allocated-data) |
| `The enclosing parent function ... for an extended __host__ __device__ lambda cannot have private or protected access within its class` | Parallel loops cannot be placed in [private or protected methods](Common_compilation_problems.md#the-enclosing-parent-function-for-an-extended-__host__-__device__-lambda-cannot-have-private-or-protected-access-within-its-class) (including GoogleTest `TEST` blocks) or in [constructors](Common_compilation_problems.md#the-enclosing-parent-function-for-an-extended-__host__-__device__-lambda-must-allow-its-address-to-be-taken) |
| `a nonstatic member reference must be relative to a specific object` | [Class members cannot be used in static functions](Common_compilation_problems.md#a-nonstatic-member-reference-must-be-relative-to-a-specific-object) |
| `X is not defined` | [Check your includes](Common_compilation_problems.md#x-is-not-defined) |
| `incomplete type is not allowed` (mentioning `IS_COVARIANT`, `IS_CONTRAVARIANT`, `Dual` or `PERIODIC`) | [Dimensions are missing required attributes](Common_compilation_problems.md#incomplete-type-is-not-allowed) |

## Problems when running the code

| Symptom | Solution |
| --- | --- |
| The code crashes or an assertion fails | Start with a [CPU debug build](../how_to/Debugging_workflow.md#cpu-debug-build) |
| My debug build does not seem to be in debug mode | [The toolchain overrides the build type](../how_to/Debugging_workflow.md#creating-a-debug-build) |
| Segmentation fault | [Look for memory misuse](../how_to/Debugging_workflow.md#segmentation-fault-or-other-memory-issues) |
| The code works on CPU but crashes on GPU | [Check where objects are accessible from](../how_to/Debugging_workflow.md#gpu-debug-build) |
| The results are different on CPU and GPU | Check that you are not [missing a synchronisation](../core_concepts/DDC_in_gyselalibxx.md#synchronicity) and that the compiler did not warn about an [implicit capture of 'this'](Common_compilation_problems.md#implicit-capture-of-this-in-extended-lambda-expression) |
| The code is slow | [Measure its performance](../how_to/profiling.md) |

## Problems with git

| Symptom | Solution |
| --- | --- |
| The submodules were not cloned | [Initialise the submodules](../how_to/Using_git.md#q-i-cloned-the-repository-but-the-submodules-were-not-cloned) |
| Git reports changes in a submodule that I did not make | [Update the submodules](../how_to/Using_git.md#q-git-reports-changes-in-the-submodule-but-i-didnt-change-this-code) |
| I committed the wrong version of a submodule | [Restore the submodule version](../how_to/Using_git.md#q-i-accidentally-committed-a-different-version-of-a-submodule-how-do-i-get-back-to-the-version-in-the-devel-branch) |

## Frequently asked questions

- [Should I create a branch on GitHub or in a private repository?](developer_FAQ.md#should-i-create-a-branch-in-the-library-on-github-or-in-a-private-repository-eg-on-gitlab)
- [What is the difference between Debug and Release mode?](developer_FAQ.md#what-is-the-difference-between-debug-and-release-mode)
- [Should I use abort, assert, or static\_assert to raise an error?](developer_FAQ.md#should-i-use-abort-assert-or-static_assert-to-raise-an-error)
