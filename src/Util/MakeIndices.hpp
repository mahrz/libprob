#ifndef _MAKE_INDICES_H_
#define _MAKE_INDICES_H_

namespace prob
{
  namespace util
  {
    /** Compile time integer arithmetic */
    namespace compile_time_list
    {
      /** Create an arbitrary integer list as type */
      template<size_t ... n>
      struct integer_list
      {
        static constexpr size_t size = sizeof...(n);

        template<size_t m>
        struct push_back
        {
          typedef integer_list<n..., m> type;
        };

        template<size_t m>
        struct push_front
        {
          typedef integer_list<m, n...> type;
        };
      };

      /** @cond PRIVATE */
      template<typename L1, typename L2>
      struct join_lists;
      /** @endcond */

      /** Join two integer_list types */
      template<template<size_t...> class L1,
      template <size_t ...> class L2,
      size_t ...n>
      struct join_lists<L1<n... >, L2<>>
      {
        typedef L1<n...> type;
      };

      /** Join two integer_list types */
      template<template<size_t...> class L1,
      template <size_t ...> class L2,
      size_t ...n,
      size_t ...m,
      size_t h>
      struct join_lists<L1<n...>, L2<h,m...>>
      {
        typedef typename join_lists<L1<n..., h>, L2<m... >>::type type;
      };

      /** Create the typed integer_list 1, ..., max */
      template<size_t max>
      struct iota_1
      {
        typedef typename iota_1<max - 1>::type::template push_back<max>::type type;
      };

      /** Create the typed integer_list 1, ..., max */
      template<>
      struct iota_1<0>
      {
        typedef integer_list<> type;
      };

      /** Create the typed integer_list 0, ..., max-1 */
      template<size_t max>
      struct iota_0
      {
        typedef typename iota_0<max - 1>::type::template push_back<max - 1>::type type;
      };

      /** Create the typed integer_list 0, ..., max-1 */
      template<>
      struct iota_0<1>
      {
        typedef integer_list<0> type;
      };

      /** Create the typed integer_list 0, ..., max-1 */
      template<>
      struct iota_0<0>
      {
        typedef integer_list<> type;
      };

      /** Create the typed integer_list n, ..., max-1 */
      template<size_t n, size_t max>
      struct iota_n
      {
        typedef typename iota_n<n, max - 1>::type::template push_back<max - 1>::type type;
      };

      /** Create the typed integer_list n, ..., max-1 */
      template<size_t n>
      struct iota_n<n, n>
      {
        typedef integer_list<> type;
      };

      /** Create the typed integer_list n, ..., max-1 */
      template<size_t n>
      struct iota_n<n, 0>
      {
        typedef integer_list<> type;
      };
    }
  }
}

#endif /* _MAKE_INDICES_H_ */
