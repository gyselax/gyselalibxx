// the implemented class (last)
#include "finitediffheatsolver.hpp"

FinitediffHeatSolver::FinitediffHeatSolver( const Configuration& config )
	: m_delta_t( config.delta_t() )
{
}

void FinitediffHeatSolver::iter( const Distributed2DField& cur, Distributed2DField& next ) const
{
	Distributed2DField::View cur_noghost = cur.noghost_view();
	Distributed2DField::View next_noghost = next.noghost_view();
	for ( int yy = 0; yy < cur_noghost.extent( DY ); ++yy ) {
		for ( int xx = 0; xx < cur_noghost.extent( DX ); ++xx ) {
			next_noghost( yy, xx ) =
					cur_noghost( yy, xx )
					+ m_delta_t* (
							(
									cur_noghost(   yy - 1, xx )
									- cur_noghost( yy,     xx ) * 2
									+ cur_noghost( yy + 1, xx )
							) / cur.delta( DY )
							+ (
									cur_noghost(   yy, xx - 1 )
									- cur_noghost( yy, xx )     * 2
									+ cur_noghost( yy, xx + 1 )
							) / cur.delta( DX )
					);
		}
	}

	// update the ghost content
	next.sync_ghosts();

	// set the time value
	next.time( cur.time() + m_delta_t );
}
