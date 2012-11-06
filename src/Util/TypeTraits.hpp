#ifndef _TYPETRAITS_H_
#define _TYPETRAITS_H_

namespace prob
{
  namespace util
  {
    /**
     * Type traits to do type level arithmetic (mainly combining lists)
     * @todo Check whether the T parameter could be removed in some of the structs
     */
    namespace traits
    {

      /** @cond PRIVATE */
      template<template<typename ...> class T, typename A, size_t n>
      struct drop;
      /** @endcond */

      /**
       * Drop the first n types of a type list
       * @tparam n The number of types to drop from the type list
       * @tparam T Type parametrising the type lists
       * @tparam A... The start of the type list
       * @tparam B The last element of the type list
       */
      template<template<typename ...> class T, typename... A, typename B, size_t n>
      struct drop<T,T<A..., B>,n>
      {
        typedef typename drop<T,T<A...>,n-1>::type type;
      };

      /**
       * Drop the first n types of a type list
       * @tparam n The number of types to drop from the type list
       * @tparam T Type parametrising the type lists
       * @tparam A... The type list
       */
      template<template<typename ...> class T, typename... A>
      struct drop<T,T<A...>,0>
      {
        typedef T<A...> type;
      };

      /** @cond PRIVATE */
      template<template<typename ...> class T, typename A, typename B>
      struct join;
      /** @endcond */

      /** 
       * Join two type lists parameterized in another type
       * This is the specialization for joining a list with void
       *
       * @tparam T Type parametrizing the type lists
       * @tparam A... The first type list
       */
      template<template<typename ...> class T, typename... A>
      struct join<T,T<A...>,void>
      {
        typedef T<A...> type;
      };

      /** 
       * Join two type lists parameterized in another type
       *
       * @tparam T Type parametrizing the type lists
       * @tparam A... The first type list
       * @tparam B... The second type list
       */
      template<template<typename ...> class T, typename... A, typename... B>
      struct join<T,T<A...>,T<B...>>
      {
        typedef T<A...,B...> type;
      };

      /** @cond PRIVATE */
      template<typename L, typename ...S>
      struct are_equivalent;
      /** @endcond */

      /** 
       * Match two type lists check whether they are equivalent
       *
       * @tparam L a type list packed into a template i.e L<T...>
       * @tparam A the head of the packed list
       * @tparam T the head of the unpacked list
       * @tparam ...B the tail of the packed list
       * @tparam ...U the tail of the unpacked list
       */
      template<typename A,
      typename T,
      typename ...B,
      typename ...U,
      template<typename ...> class L>
      struct are_equivalent<L<A,B...>, T, U...>
      {
        /**
         * True if the type lists are equivalent, false otherwise
         */
        static const bool value =
        std::is_same<typename std::decay<A>::type,T>::value &&
        are_equivalent<L<B...>,U...>::value;
      };

      /* All other cases handled from here */

      template<typename A, typename T, typename ...B,
			template<typename ...> class L>
			struct are_equivalent<L<A,B...>, T>
			{
				static const bool value =
				std::is_same<typename std::decay<A>::type,T>::value;
			};

      template<typename A, typename T, typename ... U,
			template<typename ...> class L>
			struct are_equivalent<L<A>, T, U...>
			{
				static const bool value =
				std::is_same<typename std::decay<A>::type,T>::value;
			};

      template<typename A, typename T, template<typename ...> class L>
      struct are_equivalent<L<A>, T>
      {
        static const bool value =
        std::is_same<typename std::decay<A>::type,T>::value;
      };

      /** @cond PRIVATE */
      template<size_t I, typename ...T>
      struct type_select;
      /** @endcond */

      /**
       * @brief Selects the I-th type of a type list
       */
      template<size_t I, typename HeadT, typename ...T>
      struct type_select<I, HeadT, T...>
      {
        typedef typename type_select<I-1, T...>::type type;
      };

      /**
       * @brief Selects the I-th type of a type list
       */
      template<typename HeadT, typename ...T>
      struct type_select<0, HeadT, T...>
      {
        typedef HeadT type;
      };

      /**
       * @brief Selects the I-th type of a type list
       */
      template<size_t I, typename HeadT>
      struct type_select<I, HeadT>
      {
        typedef void type;
        static_assert(true, "Index out of bounds");
      };

      /**
       * @brief Selects the I-th type of a type list
       */
      template<typename HeadT>
      struct type_select<0, HeadT>
      {
        typedef HeadT type;
      };
    }
  }
}

#endif /* _TYPETRAITS_H_ */
