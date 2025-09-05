#include <string>

#include "bumpontailequilibrium.hpp"
#include "iequilibrium.hpp"
#include "maxwellianequilibrium.hpp"

std::unique_ptr<IEquilibrium> IEquilibrium::init_from_input(
        IdxRangeSp idx_range_kinsp,
        PC_tree_t const& yaml_input_file)
{
    std::string const init_method = PCpp_string(yaml_input_file, ".Algorithm.equilibrium");

    if (init_method == "maxwellian") {
        return std::make_unique<MaxwellianEquilibrium>(
                MaxwellianEquilibrium::init_from_input(idx_range_kinsp, yaml_input_file));
    } else if (init_method == "bump_on_tail") {
        return std::make_unique<BumpontailEquilibrium>(
                BumpontailEquilibrium::init_from_input(idx_range_kinsp, yaml_input_file));
    } else {
        throw std::runtime_error("Unrecognised equilibrium requested : " + init_method);
    }
}
