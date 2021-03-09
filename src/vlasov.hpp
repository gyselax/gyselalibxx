#pragma once

#include "advection1d.hpp"
#include "distributedfield.hpp"

class Vlasov {
public:
    void operator()(const Field2D& cur, Field2D& next) const;
	
	const Advection1D& m_advec_x;
	
	const Advection1D& m_advec_vx;
};
