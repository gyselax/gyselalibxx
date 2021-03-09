#include <cassert>

// the implemented class (last)
#include "distributedfield.hpp"

using std::array;
using std::move;
using std::pair;
using std::tie;

using std::experimental::all;
using std::experimental::dynamic_extent;
using std::experimental::mdspan;
using std::experimental::subspan;

template <int N>
FieldND<N>::FieldND()
{
}

template <int N>
void FieldND<N>::swap(FieldND& other)
{
    FieldND<N> tmp = std::move(other);
    other = std::move(*this);
    *this = std::move(tmp);
}

template class FieldND<1>;
template class FieldND<2>;
