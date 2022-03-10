
#include <experimental/mdarray>

#include <iostream>
#include <iomanip>
#include <memory>

namespace stdex = std::experimental;

//================================================================================

template <
  class T,
  class ExtsA, class LayA, class AccA,
  class ExtsB, class LayB, class AccB
>
T dot_product(
  stdex::basic_mdarray<T, ExtsA, LayA, AccA> const& a,
  stdex::basic_mdarray<T, ExtsB, LayB, AccB> const& b
) //requires ExtsA::rank() == ExtsB::rank() && ExtsA::rank() == 2
{
  T result = 0;
  for(int i = 0; i < a.extent(0); ++i) {
    for(int j = 0; j < a.extent(1); ++j) {
      result += a(i, j) * b(i, j);
    }
  }
  return result;
}

//================================================================================

template <
  class T,
  class ExtsA, class LayA, class AccA
>
void fill_in_order(
  stdex::basic_mdarray<T, ExtsA, LayA, AccA>& a
) // requires ExtsA::rank() == 2
{
  T count = 0;
  for(int i = 0; i < a.extent(0); ++i) {
    for(int j = 0; j < a.extent(1); ++j) {
      a(i, j) = count++;
    }
  }
}

//================================================================================

constexpr int rows = 3;
constexpr int cols = 3;

//================================================================================

int main() {
  {
    using array_2d_dynamic = stdex::basic_mdarray<int, stdex::extents<stdex::dynamic_extent, stdex::dynamic_extent>, stdex::layout_right>;
    using array_2d_dynamic_left = stdex::basic_mdarray<int, stdex::extents<stdex::dynamic_extent, stdex::dynamic_extent>, stdex::layout_left>;

    auto a = array_2d_dynamic(rows, cols);
    auto b = array_2d_dynamic_left(rows, cols);

    fill_in_order(a);
    fill_in_order(b);

    std::cout << dot_product(a, b) << std::endl;
  }

  {
    using array_2d_10_10 = stdex::basic_mdarray<int, stdex::extents<rows, cols>, stdex::layout_right>;
    using array_2d_10_10_left = stdex::basic_mdarray<int, stdex::extents<rows, cols>, stdex::layout_right>;

    auto a = array_2d_10_10();
    auto b = array_2d_10_10_left();
    fill_in_order(a);
    fill_in_order(b);

    std::cout << dot_product(a, b) << std::endl;
  }

}