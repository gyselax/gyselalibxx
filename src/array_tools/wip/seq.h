#pragma once

#include <memory>
#include <vector>

#include "view.h"

template <class ElementType>
using SeqSpan = View1D<ElementType>;

template <class ElementType>
class Seq
{
public:
    using value_type = ElementType;

    using size_type = std::size_t;

    using difference_type = std::ptrdiff_t;

    using reference = value_type&;

    using const_reference = const value_type&;

    using pointer = value_type*;

    using const_pointer = const value_type*;

    using iterator = value_type*;

    using const_iterator = const value_type*;

    using reverse_iterator = std::reverse_iterator<iterator>;

    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

private:
    size_type m_size;

    std::unique_ptr<ElementType[]> m_data;

public:
    constexpr inline explicit Seq(size_type count)
        : m_size(count)
        , m_data(std::make_unique<ElementType>(m_size))
    {
    }

    constexpr inline Seq(size_type count, const ElementType& value)
        : m_size(count)
        , m_data(std::make_unique<ElementType>(m_size))
    {
        for (auto&& e : *this) {
            e = value;
        }
    }

    template <class ForwardIterator>
    constexpr inline Seq(ForwardIterator first, ForwardIterator last)
        : m_size(std::distance(first, last))
        , m_data(std::make_unique<ElementType>(m_size))
    {
        range_initialize(first, last);
    }

    constexpr inline Seq(const Seq& other)
        : m_size(other.m_size)
        , m_data(std::make_unique<ElementType>(m_size))
    {
        range_initialize(other.begin(), other.end());
    }

    constexpr inline Seq(Seq&& other) noexcept
        : m_size(other.m_size)
        , m_data(std::move(other.m_data))
    {
        other.m_size = 0;
        other.m_data = nullptr;
    }

    constexpr inline Seq(std::initializer_list<ElementType> init)
        : m_size(init.size())
        , m_data(std::make_unique<ElementType>(m_size))
    {
        range_initialize(init.begin(), init.end());
    }

    constexpr inline Seq& operator=(const Seq& other) noexcept {}

    constexpr inline Seq& operator=(Seq&& other) noexcept;

    constexpr inline Seq& operator=(std::initializer_list<ElementType> ilist);

    constexpr inline void assign(size_type count, const ElementType& value);

    template <class ForwardIterator>
    constexpr inline void assign(ForwardIterator first, ForwardIterator last);

    constexpr inline void assign(std::initializer_list<ElementType> ilist);

    constexpr inline reference operator[](size_type pos);

    constexpr inline const_reference operator[](size_type pos) const;

    constexpr inline reference front();

    constexpr inline const_reference front() const;

    constexpr inline reference back();

    constexpr inline const_reference back() const;

    constexpr inline ElementType* data() noexcept;

    constexpr inline const ElementType* data() const noexcept;

    constexpr inline iterator begin() noexcept;

    constexpr inline const_iterator begin() const noexcept;

    constexpr inline const_iterator cbegin() const noexcept;

    constexpr inline iterator end() noexcept;

    constexpr inline const_iterator end() const noexcept;

    constexpr inline const_iterator cend() const noexcept;

    constexpr inline reverse_iterator rbegin() noexcept;

    constexpr inline const_reverse_iterator rbegin() const noexcept;

    constexpr inline const_reverse_iterator crbegin() const noexcept;

    constexpr inline reverse_iterator rend() noexcept;

    constexpr inline const_reverse_iterator rend() const noexcept;

    constexpr inline const_reverse_iterator crend() const noexcept;

    [[nodiscard]] constexpr inline bool empty() const noexcept
    {
        return m_size == 0;
    }

    constexpr inline size_type size() const noexcept
    {
        return m_size;
    }

    constexpr inline size_type max_size() const noexcept
    {
        return m_size;
    }

    constexpr inline size_type capacity() const noexcept
    {
        return m_size;
    }

    constexpr inline void shrink_to_fit() {}

    constexpr void swap(Seq& other) noexcept
    {
        std::swap(m_data, other.m_data);
        std::swap(m_size, other.m_size);
    }

private:
    template <bool _TrivialValueTypes>
    struct UnitializedCopier
    {
        template <typename ForwardIterator>
        static void apply(ForwardIterator first, ForwardIterator last, ForwardIterator to)
        {
            for (ForwardIterator __cur = to; first != last; ++first, (void)++__cur) {
                ::new (static_cast<void*>(std::addressof(*__cur))) ElementType(*first);
            }
        }
    };

    template <>
    struct UnitializedCopier<true>
    {
        template <typename ForwardIterator>
        static void apply(ForwardIterator first, ForwardIterator last, ForwardIterator to)
        {
            return std::copy(first, last, to);
        }
    };

    template <typename ForwardIterator>
    void range_initialize(ForwardIterator first, ForwardIterator last)
    {
        typedef typename std::iterator_traits<ForwardIterator>::value_type FromElementType;
        typedef typename std::iterator_traits<ForwardIterator>::reference FromReference;
        UnitializedCopier<
                std::is_trivial_v<
                        ElementType> && std::is_trivial_v<FromElementType> && std::is_assignable_v<reference, FromReference>>::
                apply(first, last, m_data.get());
    }
};

template <class ElementType>
constexpr void swap(Seq<ElementType>& x, Seq<ElementType>& y) noexcept
{
    x.swap(y);
}
