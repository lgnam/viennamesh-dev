/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman
    Copyright (c) 2006 Dan Marsden

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_ADVANCE_IMPL_20060222_2150)
#define BOOST_FUSION_ADVANCE_IMPL_20060222_2150

namespace boost { namespace fusion {

    struct example_struct_iterator_tag;

    template<typename Struct, int Pos>
    struct example_struct_iterator;

    namespace extension
    {
        template<typename Tag>
        struct advance_impl;

        template<>
        struct advance_impl<example_struct_iterator_tag>
        {
            template<typename Iterator, typename N>
            struct apply
            {
                typedef typename Iterator::struct_type struct_type;
                typedef typename Iterator::index index;
                typedef example_struct_iterator<
                    struct_type, index::value + N::value> type;

                static type
                call(Iterator const& it)
                {
                    return type(it.struct_);
                }
            };
        };
    }
}}

#endif
