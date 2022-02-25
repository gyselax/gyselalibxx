// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <type_traits>

#include <pdi.h>

#include "ddc/chunk_span.hpp"

template <class T>
static constexpr PDI_inout_t default_access_v
        = (std::is_lvalue_reference_v<T> && !std::is_const_v<std::remove_reference_t<T>>)
                  ? PDI_INOUT
                  : PDI_OUT;

template <class T>
static constexpr PDI_inout_t chunk_default_access_v = is_writable_chunk_v<T> ? PDI_INOUT : PDI_OUT;

class PdiEvent
{
    std::string m_event;

    std::vector<std::string> m_names;

public:
    PdiEvent(std::string const& event_name) : m_event(event_name) {}

    /// @{
    /// API with access argument

    template <
            PDI_inout_t access,
            class BorrowedChunk,
            std::enable_if_t<is_borrowed_chunk_v<BorrowedChunk>, int> = 0>
    PdiEvent& with(std::string const& name, BorrowedChunk&& data)
    {
        static_assert(
                !(access & PDI_IN) || (chunk_default_access_v<BorrowedChunk> & PDI_IN),
                "Invalid access for constant data");
        auto extents = data.domain().extents().array();
        size_t rank = extents.size();
        PDI_share((name + "_rank").c_str(), &rank, PDI_OUT);
        m_names.push_back(name + "_rank");
        PDI_share((name + "_extents").c_str(), extents.data(), PDI_OUT);
        m_names.push_back(name + "_extents");
        PDI_share(name.c_str(), const_cast<chunk_value_t<BorrowedChunk>*>(data.data()), access);
        m_names.push_back(name);
        return *this;
    }

    template <
            PDI_inout_t access,
            class Arithmetic,
            std::enable_if_t<std::is_arithmetic_v<Arithmetic>, int> = 0>
    PdiEvent& with(std::string const& name, Arithmetic& data)
    {
        static_assert(
                !(access & PDI_IN) || (default_access_v<Arithmetic> & PDI_IN),
                "Invalid access for constant data");
        using value_type = std::remove_cv_t<Arithmetic>;
        PDI_share(name.c_str(), const_cast<value_type*>(&data), access);
        m_names.push_back(name);
        return *this;
    }

    template <PDI_inout_t access, class T>
    PdiEvent& and_with(std::string const& name, T&& t)
    {
        return with<access>(name, std::forward<T>(t));
    }

    /// @}
    /// API with access deduction
    /// @{

    /// Borrowed chunk overload (Chunk (const)& or ChunkSpan&& or ChunkSpan (const)&)
    template <class BorrowedChunk, std::enable_if_t<is_borrowed_chunk_v<BorrowedChunk>, int> = 0>
    PdiEvent& with(std::string const& name, BorrowedChunk&& data)
    {
        return with<chunk_default_access_v<BorrowedChunk>>(name, data);
    }

    /// Arithmetic overload (only lvalue-ref)
    template <class Arithmetic, std::enable_if_t<std::is_arithmetic_v<Arithmetic>, int> = 0>
    PdiEvent& with(std::string const& name, Arithmetic& data)
    {
        return with<default_access_v<Arithmetic>>(name, data);
    }

    /// With synonym
    template <class T>
    PdiEvent& and_with(std::string const& name, T&& t)
    {
        return with(name, std::forward<T>(t));
    }

    /// @}

    ~PdiEvent()
    {
        PDI_event(m_event.c_str());
        for (std::string const& one_name : m_names) {
            PDI_reclaim(one_name.c_str());
        }
    }
};

template <PDI_inout_t access, class DataType>
void expose_to_pdi(std::string const& name, DataType&& data)
{
    PdiEvent(name).with<access>(name, std::forward<DataType>(data));
}

template <class DataType>
void expose_to_pdi(std::string const& name, DataType&& data)
{
    PdiEvent(name).with(name, std::forward<DataType>(data));
}
