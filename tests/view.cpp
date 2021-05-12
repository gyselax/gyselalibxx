#include <gtest/gtest.h>

#include "view.h"

using namespace std;
using namespace std::experimental;

TEST(View1DTest, Constructor)
{
    std::array<double, 10> x = {0};
    std::array<double const, 10> cx = {0};
    std::array<double, 10> const ccx = {0};
    std::array<double const, 10> const ccx2 = {0};
    View1D<double> xv(x.data(), x.size());
    View1D<double> xv_(xv);
    View1D<double const> xcv(xv);
    View1D<double const> xcv_(xcv);
    View1D<double const> cxcv(cx.data(), cx.size());
    View1D<double const> cxcv_(cxcv);
    View1D<double const> ccxcv(ccx.data(), ccx.size());
    View1D<double const> ccxcv_(cxcv);
    View1D<double const> ccx2cv(ccx2.data(), ccx2.size());
    View1D<double const> ccx2cv_(ccx2cv);
}
