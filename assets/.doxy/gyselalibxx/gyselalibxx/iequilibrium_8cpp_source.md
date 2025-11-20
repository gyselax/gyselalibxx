

# File iequilibrium.cpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**initialisation**](dir_cdb336346544d0d5f695f9cdfe73a70e.md) **>** [**iequilibrium.cpp**](iequilibrium_8cpp.md)

[Go to the documentation of this file](iequilibrium_8cpp.md)


```C++
#include <string>

#include "bumpontailequilibrium.hpp"
#include "iequilibrium.hpp"
#include "maxwellianequilibrium.hpp"

std::unique_ptr<IEquilibrium> equilibrium::init_from_input(
        IdxRangeSp idx_range_kinsp,
        PC_tree_t const& yaml_input_file)
{
    std::string const init_method = PCpp_string(yaml_input_file, ".Algorithm.equilibrium");

    if (init_method == "maxwellian") {
        return std::make_unique<MaxwellianEquilibrium>(
                maxwellian_equilibrium::init_from_input(idx_range_kinsp, yaml_input_file));
    } else if (init_method == "bump_on_tail") {
        return std::make_unique<BumpontailEquilibrium>(
                bumpontail_equilibrium::init_from_input(idx_range_kinsp, yaml_input_file));
    } else {
        throw std::runtime_error("Unrecognised equilibrium requested : " + init_method);
    }
}
```


