#pragma once

// local headers
#include "distributed2dfield.hpp"

/** An interface providing access to the user-defined configuration.
 * 
 * Each configuration option is accessible through its own member function.
 */
class Configuration
{
public:
	/// The destructor
	virtual ~Configuration() = default;
	
	/** Access the number of iterations to execute
	 * @return 
	 */
	virtual int nb_iter() const = 0;
	
	/** Access the size of the global array
	 * @return the size of the global array
	 */
	virtual Coord2D global_shape() const = 0;
	
	/** Access the shape of the distribution to use
	 * @return the shape of the distribution to use
	 */
	virtual Coord2D dist_extents() const = 0;
	
	/** Access the step in time to use
	 * @return the step in time to use
	 */
	virtual double delta_t() const = 0;
	
	/** Access the step in space to use
	 * @return the step in space to use
	 */
	virtual std::array<double, 2> delta_space() const = 0;
	
};
