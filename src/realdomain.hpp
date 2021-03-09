#pragma once

#include <array>
#include <memory>

template< int N>
class RealDomainND;

template<>
class RealDomainND<1> {
	
    bool periodic;
    double xmin;
    double xmax;
    double length;
	
};

using RealDomain1D = RealDomainND<1>;

using RealDomain2D = RealDomainND<2>;

template< int N>
class RealDomainND {
	
	std::array<RealDomain1D, N> m_subgeom;
	
	std::array<std::reference_wrapper<BSplines>, N> splines() const;
	
};
