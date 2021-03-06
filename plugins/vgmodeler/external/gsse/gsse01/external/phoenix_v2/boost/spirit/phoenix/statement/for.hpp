/*=============================================================================
    Copyright (c) 2001-2004 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef PHOENIX_STATEMENT_FOR_HPP
#define PHOENIX_STATEMENT_FOR_HPP

#include <boost/spirit/phoenix/core/composite.hpp>
#include <boost/spirit/phoenix/core/compose.hpp>

namespace boost { namespace phoenix
{
    struct for_eval
    {
        template <
            typename Env
          , typename Init, typename Cond, typename Step, typename Do>
        struct apply
        {
            typedef void type;
        };

        template <
            typename RT, typename Env
          , typename Init, typename Cond, typename Step, typename Do>
        static void
        eval(
        	Env const& env
          , Init& init, Cond& cond, Step& step, Do& do_)
        {
            for (init.eval(env); cond.eval(env); step.eval(env))
                do_.eval(env);
        }
    };

    template <typename Init, typename Cond, typename Step>
    struct for_gen
    {
        for_gen(Init const& init, Cond const& cond, Step const& step)
            : init(init), cond(cond), step(step) {}

        template <typename Do>
        actor<typename as_composite<for_eval, Init, Cond, Step, Do>::type>
        operator[](Do const& do_) const
        {
            return compose<for_eval>(init, cond, step, do_);
        }

        Init init;
        Cond cond;
        Step step;
    };

    template <typename Init, typename Cond, typename Step>
    inline for_gen<Init, Cond, Step>
    for_(Init const& init, Cond const& cond, Step const& step)
    {
        return for_gen<Init, Cond, Step>(init, cond, step);
    }
}}

#endif
