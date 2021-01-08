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
#include "cartesiandistribution2d.hpp"

/** A 2D field of double distributed over a Cartesian grid.
 * 
 * A Distributed2DField represents a 2D field of double distributed in blocks
 * over a Cartesian grid described by the CartesianDistribution2D exposed by
 * the distribution() function.
 * Each local block features ghost zones that can be synchronized using the
 * sync_ghosts() function.
 * 
 * The field offer multiple "view":
 * - full_view() is a coordinate system that includes the ghosts, the point at
 *   coordinate (0,0) is in the bottom left of the ghost corner and the
 *   View::extents() includes all ghosts,
 * - noghost_view() is a coordinate system that excludes the ghosts, the point
 *   at coordinate (0,0) is in the bottom left corner of the local block and the
 *   View::extents() excludes all ghosts,
 * - ghost_view() offers a view inside each of the ghost: LEFT, RIGHT, UP and
 *   DOWN.
 * 
 * The distance between two grid points can be accessed through the delta()
 * function and the time-coordinate of the field can be set and accessed using
 * the time() function.
 */
class Distributed2DField
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
	CartesianDistribution2D m_distribution;

	/// The values
	std::vector<double> m_data;

	/// The coordinate in the time dimension for this data
	double m_time;

	/// The distance between 2 grid points in metre
	std::array<double, 2> m_delta_space;

	/// The view of the data including ghosts
	View m_full_view;

	/// The view of the data excluding ghosts
	View m_noghost_view;

	/// The views of the data in each of the ghosts
	std::array<View, 4> m_ghost_views;

	/// the MPI types to send/receive a column of ghost
	MPI_Datatype m_ghost_col;

	/// the MPI types to send/receive a row of ghost
	MPI_Datatype m_ghost_row;

public:
	/** Constructs a new Distributed2DField
	 * @param comm the MPI communicator over which the array will be distributed
	 * @param dist_shape the shape of the array distribution, i.e. the number of processes in each dimension
	 * @pre dist_shape[Dim::Y] * dist_shape[Dim::Y] must be equal to MPI_Comm_size(comm)
	 * @param global_shape the shape of the global distributed , i.e. the number of points in each dimension
	 * @pre shape[Dim::Y] must be a multiple of dist[Dim::Y]
	 * @pre shape[Dim::X] must be a multiple of dist[Dim::X]
	 * @param ghost_sizes the size of the ghosts in each dimension
	 * @param delta_space the distance between 2 grid points in metre
	 */
	Distributed2DField( MPI_Comm comm, Shape2D dist_shape, Shape2D global_shape, Shape2D ghost_sizes, std::array<double, 2> delta_space );

	/** Constructs a new Distributed2DField by copy
	 * @param other the Distributed2DField to copy
	 */
	Distributed2DField( const Distributed2DField& other ) = default;

	/** Constructs a new Distributed2DField by move
	 * @param other the Distributed2DField to move
	 */
	Distributed2DField( Distributed2DField&& other ) = default;

	/** Copy-assigns a new value to this field
	 * @param other the Distributed2DField to copy
	 * @return *this
	 */
	Distributed2DField& operator=( const Distributed2DField& other ) = default;

	/** Move-assigns a new value to this field
	 * @param other the Distributed2DField to move
	 * @return *this
	 */
	Distributed2DField& operator=( Distributed2DField&& other ) = default;

	/** Swaps this field with another
	 * @param other the Distributed2DField to swap with this one
	 */
	void swap( Distributed2DField& other );

	/** Provide a modifiable view of the data including ghosts
	 * @return a modifiable view of the data including ghosts
	 */
	View full_view() { return m_full_view; }

	/** Provide a constant view of the data including ghosts
	 * @return a constant view of the data including ghosts
	 */
	const View full_view() const { return m_full_view; }

	/** Accesses a modifiable view of the data including ghosts
	 * @param yy the coordinate in Dim Y
	 * @param xx the coordinate in Dim X
	 * @return a modifiable view of the data including ghosts
	 */
	double& full_view( int yy, int xx ) { return m_full_view( yy, xx ); }

	/** Accesses a constant view of the data including ghosts
	 * @param yy the coordinate in Dim Y
	 * @param xx the coordinate in Dim X
	 * @return a constant view of the data including ghosts
	 */
	double full_view( int yy, int xx ) const { return m_full_view( yy, xx ); }

	/** Provide a modifiable view of the data excluding ghosts
	 * @return a modifiable view of the data excluding ghosts
	 */
	View noghost_view() { return m_noghost_view; }

	/** Provide a constant view of the data excluding ghosts
	 * @return a constant view of the data excluding ghosts
	 */
	const View noghost_view() const { return m_noghost_view; }

	/** Accesses a modifiable view of the data excluding ghosts
	 * @param yy the coordinate in Dim Y
	 * @param xx the coordinate in Dim X
	 * @return a modifiable view of the data excluding ghosts
	 */
	double& noghost_view( int yy, int xx ) { return noghost_view()( yy, xx ); }

	/** Accesses a constant view of the data excluding ghosts
	 * @param yy the coordinate in Dim Y
	 * @param xx the coordinate in Dim X
	 * @return a constant view of the data excluding ghosts
	 */
	double noghost_view( int yy, int xx ) const { return noghost_view()( yy, xx ); }

	/** Provide a modifiable view of a ghost
	 * @param id the identifier of the ghost to view
	 * @return a modifiable view of a ghost
	 */
	View ghost_view( Direction2D id ) { return m_ghost_views[static_cast<int>( id )]; }

	/** Provide a constant view of a ghost
	 * @param id the identifier of the ghost to view
	 * @return a constant view of a ghost
	 */
	const View ghost_view( Direction2D id ) const { return m_ghost_views[static_cast<int>( id )]; }

	/** Accesses a modifiable view of a ghost
	 * @param id the identifier of the ghost to view
	 * @param yy the coordinate in Dim Y
	 * @param xx the coordinate in Dim X
	 * @return a modifiable view of a ghost
	 */
	double& ghost_view( Direction2D id, int yy, int xx ) { return ghost_view( id )( yy, xx ); }

	/** Accesses a constant view of a ghost
	 * @param id the identifier of the ghost to view
	 * @param yy the coordinate in Dim Y
	 * @param xx the coordinate in Dim X
	 * @return a constant view of a ghost
	 */
	double ghost_view( Direction2D id, int yy, int xx ) const { return ghost_view( id )( yy, xx ); }

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
	std::array<double, 2> delta_space() const { return m_delta_space; }

	/** Access the distance between consecutive points in metre in a given dimension
	 * @param dim the dimension to consider
	 * @return the distance between consecutive points in metre in dimension dim
	 */
	double delta( Dimension2D dim ) const { return m_delta_space[dim]; }

	/** Access the distribution of this field
	 * @return the distribution of this field
	 */
	const CartesianDistribution2D& distribution() const { return m_distribution; }

	/** Synchronize the ghosts with neighbours
	 */
	void sync_ghosts();

};
