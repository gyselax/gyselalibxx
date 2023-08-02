# Adding Documentation

There are two types of documentation:
1.  Documentation describing code structures
2.  Documentation describing general methods

## Documentation describing code structures

This documentation should be added next to the code it is describing.

Each section of code is annotated with a selection of the following:
1.  A short summary
2.  An extended summary
3.  Argument annotations
4.  Result annotations
5.  References to other relevant parts of the code

Class comments should contain a short summary and an extended summary.

Function comments should contain a short summary, an extended summary where appropriate, argument annotations, and result annotations.

Where it is helpful references to other parts of the code may also be added.

Code documentation blocks must be notated in a multi-line comment beginning with `/**`. For example:
```
/**
 * Here is a multi-line comment which will be noticed by Doxygen.
 */
```

The basic syntax is described below. For more complicated cases see the Doxygen documentation.

### Short summary

Short summaries should begin with the `@brief` tag. They should be about one sentence long.

### Extended summary

An extended summary does not need a tag. It should begin after the short summary with an empty line left between the two differentiate them.

### Argument annotation

Each argument to a function must be described. The description must begin with the `@param` tag followed by one of `[in]`, `[out]`, `[in,out]` to indicate whether the argument is constant inside the function (`[in]`), initialised by the function `[out]`, or whether the value is used by the function but the contents are then modified (`[in, out]`). This information is followed by the name of the argument and a short description separated by spaces.

### Result annotation

The result of a function must be described. The description must begin with the `@return` tag followed by a short description of what the function returns.

### References to other relevant parts of the code

If other parts of the code are relevant to your description or may be relevant in places where the code is used you can signal this with the `@see` tag which must be followed by the full name of the relevant class or function.

### Example

```cpp
/**
 * @brief A class which provides an interpolating function.
 *
 * An abstract class which implements a function allowing
 * the value of a function to be approximated at a set of
 * coordinates from a set of known values of the function.
 */
template <class DDim>
class IInterpolator
{
    using CDim = typename DDim::continuous_dimension_type;

public:
    virtual ~IInterpolator() = default;

    /**
     * @brief Approximate the value of a function at a set of coordinates using the
     * current values at a known set of interpolation points.
     *
     * @param[in, out] inout_data On input: an array containing the value of the function at the interpolation points.
     *                   On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates The coordinates where the function should be evaluated.
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     *
     * @see InterpolatorProxy::operator()
     */
    virtual ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> inout_data,
            ddc::ChunkSpan<const ddc::Coordinate<CDim>, ddc::DiscreteDomain<DDim>> coordinates)
            const = 0;
};
```

## Documentation describing general methods

Each folder containing code should also contain a `README.md` file providing an overview of the contents of the folder and any general information and references which may be helpful for understanding this code.

The `README.md` file should begin with a title:
```markdown
# <title>
```
If the folder doesn't contain a `.private` file (indicating that the folder is not yet ready to be shared publically), the new page can be referenced from any enclosing page that may exist. For example if you create a new folder in `src/`, you should add the following line to the list of folders in `src/README.md`:
```markdown
- [my_new_folder](./my_new_folder/README.md) : Short description of contents.
```

## Mathematical notation in documentation

Mathematical notation can be used in Doxygen output. It should be blocked with `$` commands. E.g:
```
$a \ne b$
```

An equation can also be printed on its own line using `$$` commands.
