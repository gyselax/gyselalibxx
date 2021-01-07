// the implemented class (last)
#include "cartesiandistribution2d.hpp"

using std::array;
using std::get;
using std::tuple;

CartesianDistribution2D::CartesianDistribution2D( MPI_Comm comm, Coord2D shape )
{
	// creation of the Cartesian communicator
	array<int, 2> cart_period = { false, false };
	MPI_Cart_create( comm, NDIM, shape.data(), cart_period.data(), true, &m_comm );
}

Coord2D CartesianDistribution2D::coord() const
{
	Coord2D result;
	Coord2D garbage;
	MPI_Cart_get( m_comm, NDIM, garbage.data(), garbage.data(), result.data() );
	return result;
}

Shape2D CartesianDistribution2D::extents() const
{
	Shape2D result;
	Coord2D garbage;
	MPI_Cart_get( m_comm, NDIM, result.data(), garbage.data(), garbage.data() );
	return result;
}

int CartesianDistribution2D::neighbour_rank( Direction2D direction )
{
	int result, garbage;
	switch ( direction ) {
	case LEFT: {
		MPI_Cart_shift( m_comm, DX, -1, &garbage, &result );
	} break;
	case RIGHT: {
		MPI_Cart_shift( m_comm, DX,  1, &garbage, &result );
	} break;
	case DOWN: {
		MPI_Cart_shift( m_comm, DY, -1, &garbage, &result );
	} break;
	case UP: {
		MPI_Cart_shift( m_comm, DY,  1, &garbage, &result );
	} break;
	}
	return result;
}
