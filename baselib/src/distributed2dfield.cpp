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

static constexpr int COMM_TAG = 1;

Distributed2DField::Distributed2DField( MPI_Comm comm, Shape2D dist_shape, Shape2D global_shape, Shape2D ghost_sizes, std::array<double, 2> delta_space )
	: m_distribution( comm, dist_shape )
	, m_data( 0 )
	, m_delta_space(delta_space)
{
	for ( int ii = 0; ii < m_data.size(); ++ii ) {
		m_data[ii] = ii;
	}
	// the shape of the local bloc without ghost
	const Shape2D noghost_shape = {
		global_shape[DY] / dist_shape[DY],
		global_shape[DX] / dist_shape[DX]
	};

	// the shape of the full local bloc, with ghost added
	const Shape2D full_shape = {
		noghost_shape[DY] + 2 * ghost_sizes[DY],
		noghost_shape[DX] + 2 * ghost_sizes[DX]
	};

	m_data.resize( full_shape[DY]*full_shape[DX] );
	mdspan<double, dynamic_extent, dynamic_extent> data2d( m_data.data(), full_shape[DY], full_shape[DX] );

	// a view of the full local bloc, including the whole data with ghosts
	m_full_view = subspan( data2d,
					pair<int, int> {0, full_shape[DY]},
					pair<int, int> {0, full_shape[DX]}
			);

	// a view of the local bloc, without ghost
	// in both dimensions, the data start after the ghost and ends before it
	m_noghost_view = subspan( data2d,
					pair<int, int> {ghost_sizes[DY], full_shape[DY] - ghost_sizes[DY]},
					pair<int, int> {ghost_sizes[DX], full_shape[DX] - ghost_sizes[DX]}
			);

	// each ghost view has the same size as the local bloc in one dimension
	// in the other dimension it contains the part of the full view excluded
	// on one of the side of m_noghost_view

	m_ghost_views[DOWN] = subspan( data2d,
					pair<int, int> {0, ghost_sizes[DY]},
					pair<int, int> {ghost_sizes[DX], full_shape[DX] - ghost_sizes[DX]}
			);

	m_ghost_views[UP] = subspan( data2d,
					pair<int, int> {ghost_sizes[DY] + noghost_shape[DY], full_shape[DY]},
					pair<int, int> {ghost_sizes[DX], full_shape[DX] - ghost_sizes[DX]}
			);

	m_ghost_views[LEFT] = subspan( data2d,
					pair<int, int> {ghost_sizes[DY], full_shape[DY] - ghost_sizes[DY]},
					pair<int, int> {0, ghost_sizes[DX]}
			);

	m_ghost_views[RIGHT] = subspan( data2d,
					pair<int, int> {ghost_sizes[DY], full_shape[DY] - ghost_sizes[DY]},
					pair<int, int> {ghost_sizes[DX] + noghost_shape[DX], full_shape[DX]}
			);

	// A type representing a column ghost (i.e. the left or right one)
	MPI_Type_vector( noghost_shape[DY], ghost_sizes[DX], full_shape[DX], MPI_DOUBLE, &m_ghost_col );
	MPI_Type_commit( &m_ghost_col );

	// A type representing a row ghost (i.e. the top or bottom one)
	MPI_Type_vector( ghost_sizes[DY], noghost_shape[DX], full_shape[DX], MPI_DOUBLE, &m_ghost_row );
	MPI_Type_commit( &m_ghost_row );
}

void Distributed2DField::sync_ghosts()
{
	// get the rank of our neighbours in the Y dimension
	int up_rank = m_distribution.neighbour_rank( UP );
	int down_rank = m_distribution.neighbour_rank( DOWN );

	// the coordinate of the last block of the same size as the ghost in the no-ghost zone
	const int LAST_Y = noghost_view().extent( DY ) - ghost_view( UP ).extent( DY );

	// send down, receive from the top
	MPI_Sendrecv(
			// send down the first row block in the no-ghost zone
			&noghost_view( 0, 0 ),      1, m_ghost_row, down_rank, COMM_TAG,
			// receive the top ghost from the top
			&ghost_view( UP, 0, 0 ),    1, m_ghost_row, up_rank,   COMM_TAG,
			m_distribution.communicator(), MPI_STATUS_IGNORE );

	// send up, receive from the bottom
	MPI_Sendrecv(
			// send up the last row block in the no-ghost zone
			&noghost_view( LAST_Y, 0 ), 1, m_ghost_row, up_rank,   COMM_TAG,
			// receive the bottom ghost from the bottom
			&ghost_view( DOWN, 0, 0 ),  1, m_ghost_row, down_rank, COMM_TAG,
			m_distribution.communicator(), MPI_STATUS_IGNORE
	);

	// get the rank of our neighbours in the X dimension
	int left_rank = m_distribution.neighbour_rank( LEFT );
	int right_rank = m_distribution.neighbour_rank( RIGHT );

	// the coordinate of the last block of the same size as the ghost in the no-ghost zone
	const int LAST_X = noghost_view().extent( DX ) - ghost_view( LEFT ).extent( DX );

	// send left, receive from the right
	MPI_Sendrecv(
			// send left the first column block in the no-ghost zone
			&noghost_view( 0, 0 ),      1, m_ghost_col, left_rank,  COMM_TAG,
			// receive the right ghost from the right
			&ghost_view( RIGHT, 0, 0 ), 1, m_ghost_col, right_rank, COMM_TAG,
			m_distribution.communicator(), MPI_STATUS_IGNORE );

	// send right, receive from the left
	MPI_Sendrecv(
			// send right the last column block in the no-ghost zone
			&noghost_view( 0, LAST_X ), 1, m_ghost_col, right_rank, COMM_TAG,
			// receive the left ghost from the left
			&ghost_view( LEFT, 0, 0 ),  1, m_ghost_col, left_rank,  COMM_TAG,
			m_distribution.communicator(), MPI_STATUS_IGNORE );
}

void Distributed2DField::swap( Distributed2DField& other )
{
	Distributed2DField tmp = std::move( other );
	other = std::move( *this );
	*this = std::move( tmp );
}
