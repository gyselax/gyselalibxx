

# File gauss\_legendre\_integration.cpp

[**File List**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**gauss\_legendre\_integration.cpp**](gauss__legendre__integration_8cpp.md)

[Go to the documentation of this file](gauss__legendre__integration_8cpp.md)


```C++
#include <array>

#include "gauss_legendre_integration.hpp"

namespace {

template <std::size_t NPoints>
constexpr std::array<long double, NPoints> s_pos;

template <std::size_t NPoints>
constexpr std::array<long double, NPoints> s_weight;

template <>
std::array<long double, 1> s_weight<1> = {2.0000000000000000000000000L};

template <>
std::array<long double, 1> s_pos<1> = {0.0000000000000000000000000L};

template <>
std::array<long double, 2>
        s_weight<2> = {1.0000000000000000000000000L, 1.0000000000000000000000000L};

template <>
std::array<long double, 2>
        s_pos<2> = {-0.5773502691896257645091488L, +0.5773502691896257645091488L};

template <>
std::array<long double, 3> s_weight<3> = {
        0.5555555555555555555555556L,
        0.8888888888888888888888889L,
        0.5555555555555555555555556L};

template <>
std::array<long double, 3> s_pos<3> = {
        -0.7745966692414833770358531L,
        +0.0000000000000000000000000L,
        +0.7745966692414833770358531L};

template <>
std::array<long double, 4> s_weight<4> = {
        0.3478548451374538573730639L,
        0.6521451548625461426269361L,
        0.6521451548625461426269361L,
        0.3478548451374538573730639L};

template <>
std::array<long double, 4> s_pos<4> = {
        -0.8611363115940525752239465L,
        -0.3399810435848562648026658L,
        +0.3399810435848562648026658L,
        +0.8611363115940525752239465L};

template <>
std::array<long double, 5> s_weight<5> = {
        0.2369268850561890875142640L,
        0.4786286704993664680412915L,
        0.5688888888888888888888889L,
        0.4786286704993664680412915L,
        0.2369268850561890875142640L};

template <>
std::array<long double, 5> s_pos<5> = {
        -0.9061798459386639927976269L,
        -0.5384693101056830910363144L,
        +0.0000000000000000000000000L,
        +0.5384693101056830910363144L,
        +0.9061798459386639927976269L};

template <>
std::array<long double, 6> s_weight<6> = {
        0.1713244923791703450402961L,
        0.3607615730481386075698335L,
        0.4679139345726910473898703L,
        0.4679139345726910473898703L,
        0.3607615730481386075698335L,
        0.1713244923791703450402961L};

template <>
std::array<long double, 6> s_pos<6> = {
        -0.9324695142031520278123016L,
        -0.6612093864662645136613996L,
        -0.2386191860831969086305017L,
        +0.2386191860831969086305017L,
        +0.6612093864662645136613996L,
        +0.9324695142031520278123016L};

template <>
std::array<long double, 7> s_weight<7> = {
        0.1294849661688696932706114L,
        0.2797053914892766679014678L,
        0.3818300505051189449503698L,
        0.4179591836734693877551020L,
        0.3818300505051189449503698L,
        0.2797053914892766679014678L,
        0.1294849661688696932706114L};

template <>
std::array<long double, 7> s_pos<7> = {
        -0.9491079123427585245261897L,
        -0.7415311855993944398638648L,
        -0.4058451513773971669066064L,
        +0.0000000000000000000000000L,
        +0.4058451513773971669066064L,
        +0.7415311855993944398638648L,
        +0.9491079123427585245261897L};

template <>
std::array<long double, 8> s_weight<8> = {
        0.1012285362903762591525314L,
        0.2223810344533744705443560L,
        0.3137066458778872873379622L,
        0.3626837833783619829651504L,
        0.3626837833783619829651504L,
        0.3137066458778872873379622L,
        0.2223810344533744705443560L,
        0.1012285362903762591525314L};

template <>
std::array<long double, 8> s_pos<8> = {
        -0.9602898564975362316835609L,
        -0.7966664774136267395915539L,
        -0.5255324099163289858177390L,
        -0.1834346424956498049394761L,
        +0.1834346424956498049394761L,
        +0.5255324099163289858177390L,
        +0.7966664774136267395915539L,
        +0.9602898564975362316835609L};

template <>
std::array<long double, 9> s_weight<9> = {
        0.0812743883615744119718922L,
        0.1806481606948574040584720L,
        0.2606106964029354623187429L,
        0.3123470770400028400686304L,
        0.3302393550012597631645251L,
        0.3123470770400028400686304L,
        0.2606106964029354623187429L,
        0.1806481606948574040584720L,
        0.0812743883615744119718922L};

template <>
std::array<long double, 9> s_pos<9> = {
        -0.9681602395076260898355762L,
        -0.8360311073266357942994298L,
        -0.6133714327005903973087020L,
        -0.3242534234038089290385380L,
        +0.0000000000000000000000000L,
        +0.3242534234038089290385380L,
        +0.6133714327005903973087020L,
        +0.8360311073266357942994298L,
        +0.9681602395076260898355762L};

template <>
std::array<long double, 10> s_weight<10> = {
        0.0666713443086881375935688L,
        0.1494513491505805931457763L,
        0.2190863625159820439955349L,
        0.2692667193099963550912269L,
        0.2955242247147528701738930L,
        0.2955242247147528701738930L,
        0.2692667193099963550912269L,
        0.2190863625159820439955349L,
        0.1494513491505805931457763L,
        0.0666713443086881375935688L};

template <>
std::array<long double, 10> s_pos<10> = {
        -0.9739065285171717200779640L,
        -0.8650633666889845107320967L,
        -0.6794095682990244062343274L,
        -0.4333953941292471907992659L,
        -0.1488743389816312108848260L,
        +0.1488743389816312108848260L,
        +0.4333953941292471907992659L,
        +0.6794095682990244062343274L,
        +0.8650633666889845107320967L,
        +0.9739065285171717200779640L};

} // namespace

template <std::size_t NPoints>
std::array<long double, NPoints> GaussLegendreCoefficients<NPoints>::pos = s_pos<NPoints>;

template <std::size_t NPoints>
std::array<long double, NPoints> GaussLegendreCoefficients<NPoints>::weight = s_weight<NPoints>;

template struct GaussLegendreCoefficients<1>;
template struct GaussLegendreCoefficients<2>;
template struct GaussLegendreCoefficients<3>;
template struct GaussLegendreCoefficients<4>;
template struct GaussLegendreCoefficients<5>;
template struct GaussLegendreCoefficients<6>;
template struct GaussLegendreCoefficients<7>;
template struct GaussLegendreCoefficients<8>;
template struct GaussLegendreCoefficients<9>;
template struct GaussLegendreCoefficients<10>;
```


