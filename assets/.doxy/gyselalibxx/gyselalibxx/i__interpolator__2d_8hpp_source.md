

# File i\_interpolator\_2d.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**i\_interpolator\_2d.hpp**](i__interpolator__2d_8hpp.md)

[Go to the documentation of this file](i__interpolator__2d_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <memory>

#include "ddc_aliases.hpp"



template <class IdxRange2D, class IdxRangeBatched>
class IInterpolator2D
{
    static_assert(IdxRange2D::rank() == 2);
    using Dim1 = typename ddc::type_seq_element_t<0, ddc::to_type_seq_t<IdxRange2D>>::
            continuous_dimension_type;
    using Dim2 = typename ddc::type_seq_element_t<1, ddc::to_type_seq_t<IdxRange2D>>::
            continuous_dimension_type;

public:
    using CoordType = Coord<Dim1, Dim2>;

    using DFieldType = DField<IdxRangeBatched>;

    using CConstFieldType = ConstField<CoordType, IdxRangeBatched>;

public:
    virtual ~IInterpolator2D() = default;

    virtual DField<IdxRangeBatched> operator()(
            DField<IdxRangeBatched> inout_data,
            ConstField<CoordType, IdxRangeBatched> coordinates) const = 0;
};



template <class IdxRange2D, class IdxRangeBatched>
class IPreallocatableInterpolator2D : public IInterpolator2D<IdxRange2D, IdxRangeBatched>
{
public:
    using typename IInterpolator2D<IdxRange2D, IdxRangeBatched>::CoordType;

public:
    virtual std::unique_ptr<IInterpolator2D<IdxRange2D, IdxRangeBatched>> preallocate() const = 0;

    DField<IdxRangeBatched> operator()(
            DField<IdxRangeBatched> inout_data,
            ConstField<CoordType, IdxRangeBatched> coordinates) const override
    {
        return (*preallocate())(inout_data, coordinates);
    }
};
```


