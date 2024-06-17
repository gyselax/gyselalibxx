#pragma once

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>
#include <species_info.hpp>

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

using CoordTor1 = ddc::Coordinate<Tor1>;
using CoordTor2 = ddc::Coordinate<Tor2>;
using CoordTor3 = ddc::Coordinate<Tor3>;
using CoordTor3D = ddc::Coordinate<Tor3, Tor2, Tor1>;

using CoordVpar = ddc::Coordinate<Vpar>;
using CoordMu = ddc::Coordinate<Mu>;
using CoordV2D = ddc::Coordinate<Mu, Vpar>;

struct GridTor1 : ddc::NonUniformPointSampling<Tor1>
{
};
struct GridTor2 : ddc::NonUniformPointSampling<Tor2>
{
};
struct GridTor3 : ddc::NonUniformPointSampling<Tor3>
{
};
struct GridVpar : ddc::NonUniformPointSampling<Vpar>
{
};
struct GridMu : ddc::NonUniformPointSampling<Mu>
{
};

using IdxTor1 = ddc::DiscreteElement<GridTor1>;
using IdxTor2 = ddc::DiscreteElement<GridTor2>;
using IdxTor3 = ddc::DiscreteElement<GridTor3>;
using IdxVpar = ddc::DiscreteElement<GridVpar>;
using IdxMu = ddc::DiscreteElement<GridMu>;

using DVecTor1 = ddc::DiscreteVector<GridTor1>;
using DVecTor2 = ddc::DiscreteVector<GridTor2>;
using DVecTor3 = ddc::DiscreteVector<GridTor3>;
using DVecVpar = ddc::DiscreteVector<GridVpar>;
using DVecMu = ddc::DiscreteVector<GridMu>;

using DDomTor1 = ddc::DiscreteDomain<GridTor1>;
using DDomTor2 = ddc::DiscreteDomain<GridTor2>;
using DDomTor3 = ddc::DiscreteDomain<GridTor3>;
using DDomVpar = ddc::DiscreteDomain<GridVpar>;
using DDomMu = ddc::DiscreteDomain<GridMu>;
using DDomTorCS = ddc::DiscreteDomain<GridTor2, GridTor1>;
using DDomTor3D = ddc::DiscreteDomain<GridTor3, GridTor2, GridTor1>;
using DDomV2D = ddc::DiscreteDomain<GridVpar, GridMu>;
using DDomSpTor3DV2D = ddc::DiscreteDomain<IDimSp, GridTor3, GridTor2, GridTor1, GridVpar, GridMu>;
using DDomSpTorCS = ddc::DiscreteDomain<IDimSp, GridTor2, GridTor1>;

template <class ElementType>
using FieldTor1 = device_t<ddc::Chunk<ElementType, DDomTor1>>;
using DFieldTor1 = FieldTor1<double>;

template <class ElementType>
using FieldTor2 = device_t<ddc::Chunk<ElementType, DDomTor2>>;
using DFieldTor2 = FieldTor2<double>;

template <class ElementType>
using FieldTor3 = device_t<ddc::Chunk<ElementType, DDomTor3>>;
using DFieldTor3 = FieldTor3<double>;

template <class ElementType>
using FieldVpar = device_t<ddc::Chunk<ElementType, DDomVpar>>;
using DFieldVpar = FieldVpar<double>;

template <class ElementType>
using FieldMu = device_t<ddc::Chunk<ElementType, DDomMu>>;
using DFieldMu = FieldMu<double>;

template <class ElementType>
using FieldTorCS = device_t<ddc::Chunk<ElementType, DDomTorCS>>;
using DFieldTorCS = FieldTorCS<double>;

template <class ElementType>
using FieldTor3D = device_t<ddc::Chunk<ElementType, DDomTor3D>>;
using DFieldTor3D = FieldTor3D<double>;

template <class ElementType>
using FieldSpTor3DV2D_host = host_t<ddc::Chunk<ElementType, DDomSpTor3DV2D>>;
using DFieldSpTor3DV2D_host = FieldSpTor3DV2D_host<double>;

template <class ElementType>
using FieldSpTor3DV2D = device_t<ddc::Chunk<ElementType, DDomSpTor3DV2D>>;
using DFieldSpTor3DV2D = FieldSpTor3DV2D<double>;

template <class ElementType>
using FieldSpTorCS_host = host_t<ddc::Chunk<ElementType, DDomSpTorCS>>;
using DFieldSpTorCS_host = FieldSpTorCS_host<double>;

template <class ElementType>
using FieldSpTorCS = device_t<ddc::Chunk<ElementType, DDomSpTorCS>>;
using DFieldSpTorCS = FieldSpTorCS<double>;

template <class ElementType>
using SpanTorCS = ddc::ChunkSpan<ElementType, DDomTorCS>;
using DSpanTorCS = SpanTorCS<double>;

template <class ElementType>
using SpanTor3D = ddc::ChunkSpan<ElementType, DDomTor3D>;
using DSpanTor3D = SpanTor3D<double>;

template <class ElementType>
using SpanSpTor3DV2D = ddc::ChunkSpan<ElementType, DDomSpTor3DV2D>;
using DSpanSpTor3DV2D = SpanSpTor3DV2D<double>;

template <class ElementType>
using ViewTor1 = ddc::ChunkView<ElementType const, DDomTor1>;
using DViewTor1 = ViewTor1<double>;

template <class ElementType>
using ViewTor2 = ddc::ChunkView<ElementType const, DDomTor2>;
using DViewTor2 = ViewTor2<double>;

template <class ElementType>
using ViewTor3 = ddc::ChunkView<ElementType const, DDomTor3>;
using DViewTor3 = ViewTor3<double>;

template <class ElementType>
using ViewVpar = ddc::ChunkView<ElementType const, DDomVpar>;
using DViewVpar = ViewVpar<double>;

template <class ElementType>
using ViewMu = ddc::ChunkView<ElementType const, DDomMu>;
using DViewMu = ViewMu<double>;

template <class ElementType>
using ViewTorCS = ddc::ChunkView<ElementType const, DDomTorCS>;
using DViewTorCS = ViewTorCS<double>;

template <class ElementType>
using ViewTor3D = ddc::ChunkView<ElementType const, DDomTor3D>;
using DViewTor3D = ViewTor3D<double>;

template <class ElementType>
using ViewSpTor3DV2D = ddc::ChunkView<ElementType const, DDomSpTor3DV2D>;
using DViewSpTor3DV2D = ViewSpTor3DV2D<double>;
