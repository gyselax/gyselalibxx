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
#include "view.hpp"

template <int N>
class FieldND {
public:
    /// ND memory view
    using View = ViewND<N>;

    /// the number of dimensions of this field
    static constexpr int NDIM = N;

private:
    /// The values
    std::vector<double> m_raw_data;

    /// The coordinate in the time dimension for this data
    double m_time;

    /// The view of the data
    View m_data;

public:
    FieldND();

    /** Constructs a new Distributed2DField by copy
   * @param other the Distributed2DField to copy
   */
    FieldND(const FieldND& other) = default;

    /** Constructs a new Distributed2DField by move
   * @param other the Distributed2DField to move
   */
    FieldND(FieldND&& other) = default;

    /** Copy-assigns a new value to this field
   * @param other the Distributed2DField to copy
   * @return *this
   */
    FieldND& operator=(const FieldND& other) = default;

    /** Move-assigns a new value to this field
   * @param other the Distributed2DField to move
   * @return *this
   */
    FieldND& operator=(FieldND&& other) = default;

    /** Swaps this field with another
   * @param other the Distributed2DField to swap with this one
   */
    void swap(FieldND& other);

    /** Provide a modifiable view of the data including ghosts
   * @return a modifiable view of the data including ghosts
   */
    View data() { return m_data; }

    /** Provide a constant view of the data including ghosts
   * @return a constant view of the data including ghosts
   */
    const View data() const { return m_data; }

    /** Sets the time for which this field is valid
   * @param new_time the time for which this field is valid
   */
    void time(double new_time) { m_time = new_time; }

    /** Access the time for which this view is valid
   * @return the time for which this view is valid
   */
    double time() const { return m_time; }
};

extern template class FieldND<1>;
using Field1D = FieldND<1>;
extern template class FieldND<2>;
using Field2D = FieldND<2>;
