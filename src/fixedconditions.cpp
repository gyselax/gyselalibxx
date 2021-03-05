// the implemented class (last)
#include "fixedconditions.hpp"

void FixedConditions::initial_condition( Data2D& data ) const
{
	// initialize everything to 0
	for ( int yy = 0; yy < data.full_view().extent( DVX ); ++yy ) {
		for ( int xx = 0; xx < data.full_view().extent( DX ); ++xx ) {
			data.full_view( yy, xx ) = 0;;
		}
	}

	// except the boundary condition on the left if our block is at the boundary itself
	if ( data.distribution().coord( DX ) == 0 ) {
		Data2D::View left_ghost = data.ghost_view( LEFT );
		for ( int yy = 0; yy < left_ghost.extent( DY ); ++yy ) {
			for ( int xx = 0; xx < data.full_view().extent( DX ); ++xx ) {
				left_ghost( yy, 0 ) = 2097.152;
			}
		}
	}

	data.time( 0 );
}
