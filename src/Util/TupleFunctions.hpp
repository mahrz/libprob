#ifndef _TUPLE_FUNCTIONS_H_
#define _TUPLE_FUNCTIONS_H_

#include "MakeIndices.hpp"

namespace prob
{
  namespace util
  {
    /** Namespace containing static functions for tuple manipulation */
    namespace tuple
    {

      /**
       * Return a subset of tuple determined by an index list
       *
       * @tparam Indices List of indices determining the subset
       * @tparam Tuple Type of the tuple
       *
       * @param idxs... Empty argument to determine the subset type
       * @param tpl Tuple input
       * @return The subset of the tuple determined by the type of idxs
       *
       */
      template<size_t ... Indices, typename Tuple>
      auto subset(const Tuple& tpl, compile_time_list::integer_list<Indices...> idxs)
      -> decltype(std::make_tuple(std::get<Indices>(tpl)...))
      {
        return std::make_tuple(std::get<Indices>(tpl)...);
      }

      /**
       * Return the tail of a tuple
       *
       * @tparam Head Type of the head element
       * @tparam Tail... Type of the tail elements
       *
       * @param tpl Tuple consisting of Head and Tail
       * @return The tail of tpl
       *
       */
      template<typename Head, typename ... Tail>
      std::tuple<Tail...> tail(const std::tuple<Head, Tail...>& tpl)
      {
        return subset(tpl, typename compile_time_list::iota_1<sizeof...(Tail)>::type());

      }

      /**
       * Return the start of a tuple
       *
       * @tparam Last Type of the last element
       * @tparam Start... Type of the start elements
       *
       * @param tpl Tuple consisting of Start and Last
       * @return The start of tpl
       *
       */
      template<typename ... T>
      auto start(const std::tuple<T...>& tpl) ->
      decltype(subset(tpl, typename compile_time_list::iota_0<sizeof...(T)-1>::type()))
      {
        return subset(tpl, typename compile_time_list::iota_0<sizeof...(T) - 1>::type());
      }

      /**
       * Return the last element of a tuple
       *
       * @tparam Head Type of the head element
       * @tparam Tail... Type of the tail elements
       *
       * @param tpl Tuple consisting of Head and Tail
       * @return The tail of tpl
       *
       */
      template<typename ... T>
      auto last(
          const std::tuple<T...>& tpl) -> decltype(std::get<sizeof...(T)-1>(tpl))
      {
        return std::get < sizeof...(T) - 1 > (tpl);
      }

      /** @cond PRIVATE */
      template<typename Indices>
      struct tuple_cons_impl;

      template<template<std::size_t...> class I,
      std::size_t... Indices
      >
      struct tuple_cons_impl<I<Indices...>>
      {
        template <
        typename A,
        typename... B,
        template <typename...> class T = std::tuple
        >
        static T<typename std::decay<A>::type, typename std::decay<B>::type...>
        tuple_cons(A && a, const T<B...>& t)
        {
          return std::make_tuple(a,std::get<Indices>(t)...);
        }

      };
      /** @endcond */

      /**
       * Constructs a new tuple from a head element and a given tail tuple
       *
       * @tparam Head Type of the head element
       * @tparam Tail... Type of the tail elements
       *
       * @param h Head element
       * @param t Tail tuple
       * @return The tuple of type std::tuple<Head, Tail...>
       *
       */
      template<typename Head, typename ... Tail,
      typename Indices = typename compile_time_list::iota_0<sizeof...(Tail)>::type,
      template<typename ...> class T = std::tuple>
      T<typename std::decay<Head>::type,
      typename std::decay<Tail>::type...>
      cons(const Head & h, const T<Tail...>& t)
      {
      	return tuple_cons_impl<Indices>::tuple_cons(h,t);
      }

      /** @cond PRIVATE */
      template<typename Indices>
      struct tuple_append_impl;

      template<template<std::size_t...> class I,
      std::size_t... Indices>
      struct tuple_append_impl<I<Indices...>>
      {
        template<
        typename A,
        typename... B,
        template <typename...> class T = std::tuple>
        static T<typename std::decay<B>::type..., typename std::decay<A>::type>
        tuple_append(A && a, const T<B...>& t)
        {
          return std::make_tuple(std::get<Indices>(t)..., a);
        }

      };
      /** @endcond */

      /**
       * Constructs a new tuple from a last element and a given tuple
       *
       * @tparam E... Type of the tuple elements
       * @tparam Last Type of the element to append
       *
       * @param tpl Tuple
       * @param l Head element
       * @return The tuple of type std::tuple<E..., Last>
       *
       */
      template<typename Last, typename ... E,
			typename Indices = typename compile_time_list::iota_0<sizeof...(E)>::type,
			template<typename ...> class T = std::tuple>
			T<typename std::decay<E>::type..., typename std::decay<Last>::type>
			append(const T<E...>& tpl, const Last & l)
			{
				return tuple_append_impl<Indices>::tuple_append(l,tpl);
			}

      /** @cond PRIVATE */
      template<typename IndicesA, typename IndicesB>
      struct tuple_concat_impl;

      template<template<std::size_t...> class I,
      std::size_t... IndicesA,
      std::size_t... IndicesB>
      struct tuple_concat_impl<I<IndicesA...>, I<IndicesB...>>
      {
        template<
        typename... A,
        typename... B,
        template <typename...> class T = std::tuple>
        static T<typename std::decay<A>::type..., typename std::decay<B>::type...>
				tuple_concat(const T<A...>& tA, const T<B...>& tB)
				{
					return std::make_tuple(std::get<IndicesA>(tA)...,std::get<IndicesB>(tB)...);
				}
      };
      /** @endcond */

      /**
       * Concats two tuples
       *
       * @tparam A... Type of the elements of the first tuple
       * @tparam B... Type of the elements of the second tuple
       *
       * @param tplA First tuple
       * @param tplB Second tuple
       * @return The tuple of type std::tuple<A..., B...>
       *
       */
      template<typename ... A, typename ... B,
			typename IndicesA = typename compile_time_list::iota_0<sizeof...(A)>::type,
			typename IndicesB = typename compile_time_list::iota_0<sizeof...(B)>::type,
			template<typename ...> class T = std::tuple>
			T<typename std::decay<A>::type..., typename std::decay<B>::type...>
			concat(const T<A...>& tplA, const T<B...>& tplB)
			{
				return tuple_concat_impl<IndicesA,IndicesB>::tuple_concat(tplA,tplB);
			}

      /** @cond PRIVATE */
      template<typename A, typename B>
      struct tuple_zip_impl;

      template<typename HeadA, typename ... TailA, typename HeadB,
          typename ... TailB>

      struct tuple_zip_impl<std::tuple<HeadA, TailA...>,
          std::tuple<HeadB, TailB...>>
      {
        static auto tuple_zip(const std::tuple<HeadA, TailA...>& tA,
            const std::tuple<HeadB, TailB...>& tB)
            -> decltype(tuple::cons(std::make_tuple(std::get<0>(tA),
                        std::get<0>(tB)),
                    tuple_zip_impl<std::tuple<TailA...>,
                    std::tuple<TailB...>>::tuple_zip(
                        tuple::tail(tA),
                        tuple::tail(tB))))
        {
          return tuple::cons(
              std::make_tuple(std::get<0> (tA), std::get<0> (tB)),
              tuple_zip_impl<std::tuple<TailA...>, std::tuple<TailB...>>::tuple_zip(
                  tuple::tail(tA), tuple::tail(tB)));
        }
      };

      template<typename HeadA, typename HeadB, typename ... TailB>
      struct tuple_zip_impl<std::tuple<HeadA>, std::tuple<HeadB, TailB...>>
      {
        static std::tuple<std::tuple<HeadA, HeadB>> tuple_zip(
            const std::tuple<HeadA>& tA, const std::tuple<HeadB, TailB...>& tB)
        {
          return std::make_tuple(
              std::make_tuple(std::get<0> (tA), std::get<0> (tB)));
        }
      };

      template<typename HeadA, typename ... TailA, typename HeadB>
      struct tuple_zip_impl<std::tuple<HeadA, TailA...>, std::tuple<HeadB>>
      {
        static std::tuple<std::tuple<HeadA, HeadB>> tuple_zip(
            const std::tuple<HeadA, TailA...>& tA, const std::tuple<HeadB>& tB)
        {
          return std::make_tuple(
              std::make_tuple(std::get<0> (tA), std::get<0> (tB)));
        }
      };

      template<typename HeadA, typename HeadB>
      struct tuple_zip_impl<std::tuple<HeadA>, std::tuple<HeadB>>
      {
        static std::tuple<std::tuple<HeadA, HeadB>> tuple_zip(
            const std::tuple<HeadA>& tA, const std::tuple<HeadB>& tB)
        {
          return std::make_tuple(
              std::make_tuple(std::get<0> (tA), std::get<0> (tB)));
        }
      };

      template<>
      struct tuple_zip_impl<std::tuple<>, std::tuple<>>
      {
        static std::tuple<std::tuple<>> tuple_zip(const std::tuple<>& tA,
            const std::tuple<>& tB)
        {
          return std::make_tuple(std::make_tuple());
        }
      };

      /** @endcond */

      /**
       * Zips two tuples to a tuple of pairs (i.e. tuples with two elements)
       *
       * @tparam A... Types of the elements of the first tuple
       * @tparam B... Types of the elements of the second tuple
       *
       * @param a First tuple
       * @param b Second tuple
       * @return The tuple of type std::tuple<std::tuple<A,B>...>
       *
       */
      template<typename ... A, typename ... B>
      auto zip(const std::tuple<A...>& a,
          const std::tuple<B...>& b)
          -> decltype(tuple_zip_impl<std::tuple<A...>,std::tuple<B...>>::tuple_zip(a,b))
      {
        return tuple_zip_impl<std::tuple<A...>, std::tuple<B...>>::tuple_zip(a,
            b);
      }

      /** @cond PRIVATE */
      template<typename F, typename A, typename ...Tail>
      struct tuple_fold_impl;

      template<typename F, typename A, typename Head, typename ... Tail>
      struct tuple_fold_impl<F, A, Head, Tail...>
      {
        static A tuple_fold(const F & f, const A & a,
            const std::tuple<Head, Tail...>& t)
        {
          return f(std::get<0> (t),
              tuple_fold_impl<F, A, Tail...>::tuple_fold(f, a, tuple::tail(t)));
        }
      };

      template<typename F, typename A, typename Head>
      struct tuple_fold_impl<F, A, Head>
      {
        static A tuple_fold(const F & f, const A & a, const std::tuple<Head>& t)
        {
          return f(std::get<0> (t), a);
        }
      };

      // Empty tuples as the tail break behaviour
      template<typename F, typename A>
      struct tuple_fold_impl<F, A, std::tuple<>>
      {
        static A tuple_fold(const F & f, const A & a,
            const std::tuple<std::tuple<>>& t)
        {
          return a;
        }
      };

      template<typename F, typename A>
      struct tuple_fold_impl<F, A>
      {
        static A tuple_fold(const F & f, const A & a, const std::tuple<>& t)
        {
          return a;
        }
      };

      /** @endcond */

      /**
       * Folds a tuple of elements of the same type (or same methods)
       *
       * @tparam F Type of the fold function (needs to be (B,A) -> A)
       * @tparam A Type of the initial and accumulation value
       * @tparam B... Type of the elements the tuple
       *
       * @param f The fold function
       * @param a The initial value
       * @param tpl The tuple
       * @return The folded value of the tuple
       *
       */
      template<typename F, typename A, typename ... B>
      A fold(const F & f, const A & a, const std::tuple<B...>& tpl)
      {
        return tuple_fold_impl<F, A, B...>::tuple_fold(f, a, tpl);
      }
    }
  }
}

#endif /* _TUPLE_FUNCTIONS_H_ */
