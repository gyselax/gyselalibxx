#pragma once

// the MPI library must be included before C++ headers
#include <mpi.h>

// local headers
#include "coord2d.hpp"

class CartesianDistribution2D
{
	/// the number of dimensions of this distribution
	static constexpr int NDIM = 2;

	/// the MPI Cartesian communicator over which the array is distributed
	MPI_Comm m_comm;

public:
	/** Constructs a new CartesianDistribution2D
	 * @param comm the communicator supporting the distribution
	 * @param shape the shape (extents) of the distribution, i.e. number of processes in each dimension
	 */
	CartesianDistribution2D( MPI_Comm comm, Shape2D shape );

	/** Accesses the shape of the distribution
	 * @return the shape of the distribution
	 */
	Shape2D extents() const;

	/** Accesses the shape of the distribution in a given dimension
	 * @param dim the dimension to consider
	 * @return the shape of the distribution in the requested dimension
	 */
	int extent( Dimension2D dim ) const { return extents()[dim]; }

	/** Accesses the local coordinate in the distribution
	 * @return the local coordinate in the distribution
	 */
	Coord2D coord() const;

	/** Accesses the local coordinate in the distribution in a given dimension
	 * @param dim the dimension to consider
	 * @return the local coordinate in the distribution in the requested dimension
	 */
	int coord( Dimension2D dim ) const { return coord()[dim]; }

	/** Accesses the MPI communicator supporting the distribution
	 * @return the MPI communicator supporting the distribution
	 */
	const MPI_Comm communicator() const { return m_comm; }

	/** Accesses the MPI communicator supporting the distribution
	 * @return the MPI communicator supporting the distribution
	 */
	MPI_Comm communicator() { return m_comm; }

	/** Accesses the rank of the neighbour process in a given direction
	 * @return the rank of the neighbour process in a given direction
	 */
	int neighbour_rank( Direction2D direction );

	/** Accesses the size of the MPI communicator supporting the distribution
	 * @return the size of the MPI communicator supporting the distribution
	 */
	int size() const { int result; MPI_Comm_size( m_comm, &result ); return result; }

	/** Accesses the local rank in the MPI communicator supporting the distribution
	 * @return the local rank in the MPI communicator supporting the distribution
	 */
	int rank() const { int result; MPI_Comm_rank( m_comm, &result ); return result; }

};
