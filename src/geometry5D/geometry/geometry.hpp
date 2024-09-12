// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "species_info.hpp"

/*
 * @file geometry.hpp
 *
 */

/**
 * @brief Define non periodic radial dimension @f$r@f$.
 */
struct Tor1
{
    /**
     * The periodicity of the radial dimension.
     * This is a compile time constant which allows the code to avoid
     * unnecessary if conditions.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define periodic poloidal dimension @f$\theta@f$.
 */
struct Tor2
{
    /**
     * The periodicity of the poloidal dimension.
     * This is a compile time constant which allows the code to avoid
     * unnecessary if conditions.
     */
    static bool constexpr PERIODIC = true;
};
/**
 * @brief Define periodic toroidal dimension @f$\varphi@f$.
 */
struct Tor3
{
    /**
     * The periodicity of the toroidal dimension.
     * This is a compile time constant which allows the code to avoid
     * unnecessary if conditions.
     */
    static bool constexpr PERIODIC = true;
};
/**
 * @brief Define non periodic parallel velocity @f$v_\parallel@f$.
 */
struct Vpar
{
    /**
     * The periodicity of the parallel velocity.
     * This is a compile time constant which allows the code to avoid
     * unnecessary if conditions.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define non periodic magnetic momentum @f$\mu@f$.
 */
struct Mu
{
    /**
     * The periodicity of the magnetic momentum.
     * This is a compile time constant which allows the code to avoid
     * unnecessary if conditions.
     */
    static bool constexpr PERIODIC = false;
};

using CoordTor1 = Coord<Tor1>;
using CoordTor2 = Coord<Tor2>;
using CoordTor3 = Coord<Tor3>;
using CoordTor3D = Coord<Tor3, Tor2, Tor1>;

using CoordVpar = Coord<Vpar>;
using CoordMu = Coord<Mu>;
using CoordV2D = Coord<Mu, Vpar>;

struct GridTor1 : NonUniformGridBase<Tor1>
{
};
struct GridTor2 : NonUniformGridBase<Tor2>
{
};
struct GridTor3 : NonUniformGridBase<Tor3>
{
};
struct GridVpar : NonUniformGridBase<Vpar>
{
};
struct GridMu : NonUniformGridBase<Mu>
{
};

using IdxTor1 = Idx<GridTor1>;
using IdxTor2 = Idx<GridTor2>;
using IdxTor3 = Idx<GridTor3>;
using IdxVpar = Idx<GridVpar>;
using IdxMu = Idx<GridMu>;

using IdxStepTor1 = IdxStep<GridTor1>;
using IdxStepTor2 = IdxStep<GridTor2>;
using IdxStepTor3 = IdxStep<GridTor3>;
using IdxStepVpar = IdxStep<GridVpar>;
using IdxStepMu = IdxStep<GridMu>;

using IdxRangeTor1 = IdxRange<GridTor1>;
using IdxRangeTor2 = IdxRange<GridTor2>;
using IdxRangeTor3 = IdxRange<GridTor3>;
using IdxRangeVpar = IdxRange<GridVpar>;
using IdxRangeMu = IdxRange<GridMu>;
using IdxRangeTorCS = IdxRange<GridTor2, GridTor1>;
using IdxRangeTor3D = IdxRange<GridTor3, GridTor2, GridTor1>;
using IdxRangeV2D = IdxRange<GridVpar, GridMu>;
using IdxRangeSpTor3DV2D = IdxRange<Species, GridTor3, GridTor2, GridTor1, GridVpar, GridMu>;
using IdxRangeSpTorCS = IdxRange<Species, GridTor2, GridTor1>;

template <class ElementType>
using FieldMemTor1 = FieldMem<ElementType, IdxRangeTor1>;
using DFieldMemTor1 = FieldMemTor1<double>;

template <class ElementType>
using FieldMemTor2 = FieldMem<ElementType, IdxRangeTor2>;
using DFieldMemTor2 = FieldMemTor2<double>;

template <class ElementType>
using FieldMemTor3 = FieldMem<ElementType, IdxRangeTor3>;
using DFieldMemTor3 = FieldMemTor3<double>;

template <class ElementType>
using FieldMemVpar = FieldMem<ElementType, IdxRangeVpar>;
using DFieldMemVpar = FieldMemVpar<double>;

template <class ElementType>
using FieldMemMu = FieldMem<ElementType, IdxRangeMu>;
using DFieldMemMu = FieldMemMu<double>;

template <class ElementType>
using FieldMemTorCS = FieldMem<ElementType, IdxRangeTorCS>;
using DFieldMemTorCS = FieldMemTorCS<double>;

template <class ElementType>
using FieldMemTor3D = FieldMem<ElementType, IdxRangeTor3D>;
using DFieldMemTor3D = FieldMemTor3D<double>;

template <class ElementType>
using FieldSpTor3DV2D_host = host_t<FieldMem<ElementType, IdxRangeSpTor3DV2D>>;
using DFieldSpTor3DV2D_host = FieldSpTor3DV2D_host<double>;

template <class ElementType>
using FieldMemSpTor3DV2D = FieldMem<ElementType, IdxRangeSpTor3DV2D>;
using DFieldMemSpTor3DV2D = FieldMemSpTor3DV2D<double>;

template <class ElementType>
using FieldSpTorCS_host = host_t<FieldMem<ElementType, IdxRangeSpTorCS>>;
using DFieldSpTorCS_host = FieldSpTorCS_host<double>;

template <class ElementType>
using FieldMemSpTorCS = FieldMem<ElementType, IdxRangeSpTorCS>;
using DFieldMemSpTorCS = FieldMemSpTorCS<double>;

template <class ElementType>
using FieldTorCS = Field<ElementType, IdxRangeTorCS>;
using DFieldTorCS = FieldTorCS<double>;

template <class ElementType>
using FieldTor3D = Field<ElementType, IdxRangeTor3D>;
using DFieldTor3D = FieldTor3D<double>;

template <class ElementType>
using FieldSpTor3DV2D = Field<ElementType, IdxRangeSpTor3DV2D>;
using DFieldSpTor3DV2D = FieldSpTor3DV2D<double>;

template <class ElementType>
using ConstFieldTor1 = ConstField<ElementType const, IdxRangeTor1>;
using DConstFieldTor1 = ConstFieldTor1<double>;

template <class ElementType>
using ConstFieldTor2 = ConstField<ElementType const, IdxRangeTor2>;
using DConstFieldTor2 = ConstFieldTor2<double>;

template <class ElementType>
using ConstFieldTor3 = ConstField<ElementType const, IdxRangeTor3>;
using DConstFieldTor3 = ConstFieldTor3<double>;

template <class ElementType>
using ConstFieldVpar = ConstField<ElementType const, IdxRangeVpar>;
using DConstFieldVpar = ConstFieldVpar<double>;

template <class ElementType>
using ConstFieldMu = ConstField<ElementType const, IdxRangeMu>;
using DConstFieldMu = ConstFieldMu<double>;

template <class ElementType>
using ConstFieldTorCS = ConstField<ElementType const, IdxRangeTorCS>;
using DConstFieldTorCS = ConstFieldTorCS<double>;

template <class ElementType>
using ConstFieldTor3D = ConstField<ElementType const, IdxRangeTor3D>;
using DConstFieldTor3D = ConstFieldTor3D<double>;

template <class ElementType>
using ConstFieldSpTor3DV2D = ConstField<ElementType const, IdxRangeSpTor3DV2D>;
using DConstFieldSpTor3DV2D = ConstFieldSpTor3DV2D<double>;
