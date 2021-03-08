#include <cassert>

// the implemented class (last)
#include "distributed2dfield.hpp"

using std::array;
using std::move;
using std::pair;
using std::tie;

using std::experimental::all;
using std::experimental::dynamic_extent;
using std::experimental::mdspan;
using std::experimental::subspan;

template<int N>
DataND<N>::DataND ( Shape<N> global_shape, MPI_Comm comm, Dimension distributed_dim, std::array<double, N> delta )
	: m_distribution( comm )
	, m_distributed_dim(distributed_dim)
	, m_raw_data( 0 )
	, m_delta_space( delta )
{
	for ( int ii = 0; ii < m_raw_data.size(); ++ii ) {
		m_raw_data[ii] = ii;
	}

	// the shape of the full local bloc, with ghost added
	int mpi_size;
	MPI_Comm_size(comm, &mpi_size);
	Shape<N> full_shape = global_shape;
	assert(full_shape[distributed_dim] % mpi_size ==0);
	full_shape[distributed_dim] /= mpi_size;
	

	m_raw_data.resize( full_shape[DVX]*full_shape[DX] );
	mdspan<double, dynamic_extent, dynamic_extent> data2d( m_raw_data.data(), full_shape[DVX], full_shape[DX] );

	// a view of the full local bloc, including the whole data with ghosts
	m_data = subspan( data2d,
					pair<int, int> {0, full_shape[DVX]},
					pair<int, int> {0, full_shape[DX]}
			);
}

template<int N>
void DataND<N>::swap( DataND& other )
{
	DataND<N> tmp = std::move( other );
	other = std::move( *this );
	*this = std::move( tmp );
}

template class DataND<1>;
template class DataND<2>;
