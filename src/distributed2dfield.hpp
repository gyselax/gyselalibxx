#pragma once

// the MPI library must be included before C++ headers
#include <mpi.h>

// standard C++ library headers
#include <array>
#include <vector>

// library headers
#include <experimental/mdspan>

// local headers
#include "coord2d.hpp"

template < int N >
class DataND
{
public:
	/// A type for a 2D view of the memory
	using View = std::experimental::basic_mdspan<double,
		  std::experimental::extents<std::experimental::dynamic_extent, std::experimental::dynamic_extent>,
		  std::experimental::layout_stride<std::experimental::dynamic_extent, 1> >;

	/// the number of dimensions of this field
	static constexpr int NDIM = 2;

private:
	/// The descriptor of data distribution over MPI
	MPI_Comm m_distribution;

	/// The descriptor of data distribution over MPI
	Dimension m_distributed_dim;

	/// The values
	std::vector<double> m_raw_data;

	/// The coordinate in the time dimension for this data
	double m_time;

	/// The distance between 2 grid points in metre
	std::array<double, N> m_delta_space;

	/// The view of the data including ghosts
	View m_data;

public:
	DataND ( Shape<N> global_shape, MPI_Comm comm, Dimension distributed_dim, std::array<double, N> delta );

	/** Constructs a new Distributed2DField by copy
	 * @param other the Distributed2DField to copy
	 */
	DataND ( const DataND& other ) = default;

	/** Constructs a new Distributed2DField by move
	 * @param other the Distributed2DField to move
	 */
	DataND ( DataND&& other ) = default;

	/** Copy-assigns a new value to this field
	 * @param other the Distributed2DField to copy
	 * @return *this
	 */
	DataND& operator=( const DataND& other ) = default;

	/** Move-assigns a new value to this field
	 * @param other the Distributed2DField to move
	 * @return *this
	 */
	DataND& operator=( DataND&& other ) = default;

	/** Swaps this field with another
	 * @param other the Distributed2DField to swap with this one
	 */
	void swap( DataND& other );

	/** Provide a modifiable view of the data including ghosts
	 * @return a modifiable view of the data including ghosts
	 */
	View data() { return m_data; }

	/** Provide a constant view of the data including ghosts
	 * @return a constant view of the data including ghosts
	 */
	const View data() const { return m_data; }

	/** Accesses a modifiable view of the data including ghosts
	 * @param yy the coordinate in Dim Y
	 * @param xx the coordinate in Dim X
	 * @return a modifiable view of the data including ghosts
	 */
	double& data( int yy, int xx ) { return m_data( yy, xx ); }

	/** Accesses a constant view of the data including ghosts
	 * @param yy the coordinate in Dim Y
	 * @param xx the coordinate in Dim X
	 * @return a constant view of the data including ghosts
	 */
	double data( int yy, int xx ) const { return m_data( yy, xx ); }

	/** Sets the time for which this field is valid
	 * @param new_time the time for which this field is valid
	 */
	void time( double new_time ) { m_time = new_time; }

	/** Access the time for which this view is valid
	 * @return the time for which this view is valid
	 */
	double time() const { return m_time; }

	/** Access the distance between consecutive points in metre
	 * @return the distance between consecutive points in metre
	 */
	std::array<double, N> delta_space() const { return m_delta_space; }

	/** Access the distance between consecutive points in metre in a given dimension
	 * @param dim the dimension to consider
	 * @return the distance between consecutive points in metre in dimension dim
	 */
	double delta( Dimension dim ) const { return m_delta_space[dim]; }

	/** Access the distribution of this field
	 * @return the distribution of this field
	 */
	MPI_Comm distribution() const { return m_distribution; }

};

extern template class DataND<1>;
extern template class DataND<2>;
