#pragma once

#include <string>
#include <type_traits>

#include <pdi.h>

#include "ddc/chunck_span.hpp"

template <class ElementType>
static constexpr PDI_inout_t default_access_v = std::is_const_v<ElementType> ? PDI_OUT : PDI_INOUT;

template <class ElementType, class SupportType, class LayoutStridedPolicy>
static constexpr PDI_inout_t default_access_v<ChunkSpan<
        ElementType,
        SupportType,
        LayoutStridedPolicy> const&> = default_access_v<ElementType>;

template <class ElementType>
static constexpr bool is_chunckspan_v = false;

template <class ElementType, class SupportType, class LayoutStridedPolicy>
static constexpr bool
        is_chunckspan_v<ChunkSpan<ElementType, SupportType, LayoutStridedPolicy>> = true;

template <class ElementType, class SupportType>
static constexpr bool is_chunckspan_v<Chunk<ElementType, SupportType>> = true;

class PdiEvent
{
    std::string m_event;

    std::vector<std::string> m_names;

public:
    PdiEvent(std::string const& event_name) : m_event(event_name) {}

    template <PDI_inout_t access, class ElementType, class SupportType, class LayoutStridedPolicy>
    PdiEvent& with(
            std::string const& name,
            ChunkSpan<ElementType, SupportType, LayoutStridedPolicy> const& data)
    {
        static_assert(
                !(access & PDI_IN) || !std::is_const_v<ElementType>,
                "Invalid access for constant data");
        auto extents = data.domain().extents().array();
        size_t rank = extents.size();
        PDI_share((name + "_rank").c_str(), &rank, PDI_OUT);
        m_names.push_back(name + "_rank");
        PDI_share(
                (name + "_extents").c_str(),
                const_cast<DiscreteVectorElement*>(extents.data()),
                PDI_OUT);
        m_names.push_back(name + "_extents");
        PDI_share(name.c_str(), const_cast<std::remove_const_t<ElementType>*>(data.data()), access);
        m_names.push_back(name);
        return *this;
    }

    template <PDI_inout_t access, class ElementType>
    std::enable_if_t<
            !is_chunckspan_v<std::remove_cv_t<std::remove_reference_t<ElementType>>>,
            PdiEvent>&
    with(std::string const& name, ElementType& data)
    {
        static_assert(
                !(access & PDI_IN) || !std::is_const_v<ElementType>,
                "Invalid access for constant data");
        PDI_share(name.c_str(), const_cast<std::remove_const_t<ElementType>*>(&data), access);
        m_names.push_back(name);
        return *this;
    }

    template <class ElementType, class SupportType, class LayoutStridedPolicy>
    PdiEvent& with(
            std::string const& name,
            ChunkSpan<ElementType, SupportType, LayoutStridedPolicy> const& data)
    {
        return with<default_access_v<ElementType>>(name, data);
    }

    template <class ElementType>
    PdiEvent& with(std::string const& name, ElementType& data)
    {
        return with<default_access_v<ElementType>>(name, data);
    }

    template <class ElementType>
    PdiEvent& and_with(std::string const& name, ElementType& data)
    {
        return with(name, data);
    }

    template <PDI_inout_t access, class ElementType>
    PdiEvent& and_with(std::string const& name, ElementType& data)
    {
        return with<default_access_v<ElementType>>(name, data);
    }

    ~PdiEvent()
    {
        PDI_event(m_event.c_str());
        for (auto&& one_name : m_names) {
            PDI_reclaim(one_name.c_str());
        }
    }
};

template <PDI_inout_t access, class ElementType>
void expose_to_pdi(std::string const& name, ElementType& data)
{
    PdiEvent(name).with<access>(name, data);
}

template <class ElementType>
void expose_to_pdi(std::string const& name, ElementType& data)
{
    PdiEvent(name).with(name, data);
}
