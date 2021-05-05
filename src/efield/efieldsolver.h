#pragma once

#include "iefieldsolver.h"

class EfieldSolver : public IEfieldSolver
{
public:
    void operator() ( DBlockViewX& ex, const DBlockViewXVx& fdistribu ) const override;
};
