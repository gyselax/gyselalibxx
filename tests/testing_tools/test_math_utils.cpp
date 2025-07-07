

KOKKOS_FUNCTION double cyl_bessel_j()
{
    double fct = 1;
    double sum = 0;
    for (int k = 0; k < 10; fct *= ++k) {
        sum += std::pow(-1, k) * std::pow(x / 2, 2 * k) / std::pow(fct, 2);
    }
    return sum;
}
