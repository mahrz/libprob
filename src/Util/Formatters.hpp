#ifndef _FORMATTERS_H_
#define _FORMATTERS_H_

#include <iostream>
#include <cstddef>
#include "../RandomVariable.hpp"
/** @cond PRIVATE */

template<std::size_t> struct int_{};

/* Pretty print tuples, space separated; used by the distribution printer */
template<class Ch, class Tr, class T, std::size_t I>
void print_tuple(std::basic_ostream<Ch,Tr>& os, T const& t, int_<I>)
{
  print_tuple(os, t, int_<I-1>());
  os << " " << std::get<I>(t);
}

template<class Ch, class Tr, class T>
void print_tuple(std::basic_ostream<Ch,Tr>& os, T const& t, int_<0>)
{
  os << std::get<0>(t);
}

template<class Ch, class Tr, class T>
void print_tuple(std::basic_ostream<Ch,Tr>& os, T const& t, int_<-1>)
{

}

template<class Ch, class Tr, class... Args>
std::ostream& operator<<(std::basic_ostream<Ch,Tr>& os,
			 std::tuple<Args...> const& t)
{
  print_tuple(os, t, int_<sizeof...(Args)-1>());
  return os;
}

/* Print random_events simply as their values */
template<class Ch, class Tr, class... Args>
std::ostream& operator<<(std::basic_ostream<Ch,Tr>& os,
			 prob::random_event const& r)
{
  return os << r._val;
}

template<class Ch, class Tr, class... Args>
std::ostream& operator<<(std::basic_ostream<Ch,Tr>& os,
       prob::given const& g)
{
  return os << "|";
}
/** @endcond */

#endif /* _FORMATTERS_H_ */
