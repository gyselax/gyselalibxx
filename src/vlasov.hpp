#pragma once

#include "advection1d.hpp"
#include "distributedfield.hpp"
#include "geometry.hpp"

class Vlasov {
public:
	const Geometry& m_geom;
	
    void operator()(const Field2D& cur, Field2D& next) const;
	
	const Advection1D& m_advec_x;
	
	const Advection1D& m_advec_vx;
	
};
