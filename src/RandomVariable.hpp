#ifndef _RANDOMVARIABLE_H_
#define _RANDOMVARIABLE_H_

#include "prob.hpp"

/**
 * @file RandomVariable.hpp
 *
 * @brief Declaration of random variables and events.
 *
 * Macros to declare random variables, classes to deal
 * with random variable events and some internals to extract static
 * extents.
 */

/**
 * @defgroup RVAR Random Variables
 *
 * @brief Declare random variable types and index random variable events.
 *
 * @{
 */

/**
 * @def RVAR Define a random variable type.
 *
 * @param x The name of the random variable.
 */
#define RVAR(x) struct x : public prob::random_event \
  { x(int ix) : prob::random_event(ix) {}		\
    x() : prob::random_event() {}	\
    x(const random_event& other) : prob::random_event(other._val) {};	\
    template<typename B> \
    prob::_given<x,B> operator |(const B& b) const	\
    { return prob::_given<x,B>(*this,b);}		\
    static constexpr int extent() { return 1; };	\
    static constexpr bool static_rvar() { return false; }; \
    static std::string label() { return std::string(#x); };	\
    };

#define RVAR_STATIC(x,e) struct x : public prob::random_event \
  { x(int ix) : prob::random_event(ix) {}	\
    x() : prob::random_event() {}	\
    x(const random_event& other) : prob::random_event(other._val) {}; \
    template<typename B> \
      prob::_given<x,B> operator |(const B& b) const \
    { return prob::_given<x,B>(*this,b);}	\
    static constexpr int extent() { return e; }; \
    static constexpr bool static_rvar() { return true; }; \
    static std::string label() { return std::string(#x); };	\
  };

/** @} */

namespace prob
{
	/**
	 * @addtogroup RVAR
	 *
	 * @{
	 */

  /**
   * Dummy type to separate conditioned variables from conditional variables.
   * Is used when lambdas are supplied
   *
   * @code
   *  p.each_index( [&] (const& A, given, const& B) { return; } );
   * @endcode
   *
   */
  struct given
  {
    /**
     * Constructor argument is ignored, but this allows to construct
     * an object of type given the same as any random variable. This
     * is needed for internal indexing.
     */
    given(int)
    {
    }

    /**
     * Empty constructor
     */
    given()
    {
    }
  };

  /** @} */

  /**
   * Alternative type to split conditioned variables from conditional
   * variables, already incorporating the last conditioned and the first
   * conditional variable. This is implicitly used when the | operator is used
   * to concatenate conditioned and conditional variables.
   */
  template<typename A, typename B>
  struct _given
  {
    _given(const A& a, const B& b) :
        _a(a), _b(b)
    {
    }

    A _a;
    B _b;
  };

  namespace core
  {
    template<typename T>
    struct extent_assert
    {
      static void assert_static()
      {
        static_assert(!T::static_rvar(),
        		"You are using a staticly sized random variable in a dynamic context");
      }

      static void assert_dynamic()
      {
        static_assert(T::static_rvar(),
        		"You are using a dynamicly sized random variable in a static context");
      }

      static bool is_static()
      {
      	return T::static_rvar();
      }

    };

    template<>
    struct extent_assert<given>
    {
      static void assert_static()
      {
      }
      static void assert_dynamic()
      {
      }

      static bool is_static()
      {
         return false;
      }
    };

    template<typename A, typename B>
    struct extent_assert<_given<A, B>>
    {
      static void assert_static()
      {
        extent_assert<A>::assert_static();
        extent_assert<B>::assert_static();
      }

      static void assert_dynamic()
      {
        extent_assert<A>::assert_dynamic();
        extent_assert<B>::assert_dynamic();
      }

      static bool is_static()
      {
         return extent_assert<A>::is_static() ||
        		 	 	extent_assert<B>::is_static();
      }
    };
  }

  /** @cond PRIVATE */
  template<typename ... T>
  struct extents_assert;
  /** @endcond */

  /**
   * Assert on a template parameter list whether all
   * types are static or dynamic random variables.
   * Based on a runtime static assertion (failure
   * will result in a compiler error!).
   */
  template<typename Head, typename ...Tail>
  struct extents_assert<Head, Tail...>
  {
    static void assert_static()
    {
      core::extent_assert<Head>::assert_static();
      extents_assert<Tail...>::assert_static();
    }

    static void assert_dynamic()
    {
      core::extent_assert<Head>::assert_dynamic();
      extents_assert<Tail...>::assert_dynamic();
    }

    static bool is_static()
		{
			 return core::extent_assert<Head>::is_static() ||
							extents_assert<Tail...>::is_static();
		}
  };

  /**
   * Assert on a template parameter whether the
   * types is a static or dynamic random variable.
   * Based on a runtime static assertion (failure
   * will result in a compiler error!).
   */
  template<typename Head>
  struct extents_assert<Head>
  {
    static void assert_static()
    {
      core::extent_assert<Head>::assert_static();
    }

    static void assert_dynamic()
    {
      core::extent_assert<Head>::assert_dynamic();
    }

    static bool is_static()
		{
			 return core::extent_assert<Head>::is_static();
		}
  };

	/**
	 * @addtogroup RVAR
	 *
	 * @{
	 */

  /**
   * In some cases wrapping an index in a random event of the corresponding
   * random variable class is a bit clumsy, therefore indices are always
   * accessed using the read_index wrapper, so in the case that the type
   * is unambiguous integers can be used instead of random_events
   */
  template<typename T>
  struct read_index
  {
    static int read(const T& i)
    {
      return i._val;
    }
    ;
    static int read(T&& i)
    { return i._val;};
  };

  /**
   * Integer specialization
   */
  template<>
  struct read_index<int>
  {
    static int read(const int& i)
    {
      return i;
    }
    ;
  };

  /**
   * Classes deriving random_event are denoting different random
   * variables. Different instantiations are then used to differentiate
   * different events indexed by integers.
   */
  struct random_event
  {
    random_event(int val) :
        _val(val)
    {
    };

    random_event(const random_event& other) :
        _val(other._val)
    {
    };

    random_event() :
        _val(0)
    {
    };

    int _val;

    template<typename B>
    _given<random_event, B> operator |(const B& b) const
    {
      return _given<random_event, B>(*this, b);
    }

    template<typename B>
    random_event operator +(const B& b) const
    {
      return random_event(_val + read_index<B>::read(b));
    }

    template<typename B>
    random_event operator -(const B& b) const
    {
      return random_event(_val - read_index<B>::read(b));
    }

    template<typename B>
    random_event operator *(const B& b) const
    {
      return random_event(_val * read_index<B>::read(b));
    }

    template<typename B>
    random_event operator /(const B& b) const
    {
      return random_event(_val / read_index<B>::read(b));
    }

    template<typename B>
    random_event operator %(const B& b) const
    {
      return random_event(_val % read_index<B>::read(b));
    }

    random_event&
    operator++()
    {
      ++_val;
      return *this;
    }

    random_event operator++(int)
    {
      random_event re = *this;
      ++_val;
      return re;
    }

    random_event&
    operator--()
    {
      --_val;
      return *this;
    }

    random_event operator--(int)
    {
      random_event re = *this;
      --_val;
      return re;
    }

    template<typename B>
		random_event operator +=(const B& b)
		{
      _val += read_index<B>::read(b);
      return *this;
		}

    template<typename B>
		random_event operator -=(const B& b)
		{
      _val -= read_index<B>::read(b);
      return *this;
		}

    template<typename B>
		random_event operator *=(const B& b)
		{
      _val *= read_index<B>::read(b);
      return *this;
		}

    template<typename B>
		random_event operator /=(const B& b)
		{
      _val /= read_index<B>::read(b);
      return *this;
		}

    template<typename B>
		random_event operator %=(const B& b)
		{
      _val %= read_index<B>::read(b);
      return *this;
		}

    template<typename B>
		bool operator ==(const B& b) const
		{
      return _val == read_index<B>::read(b);
		}

    template<typename B>
		bool operator !=(const B& b) const
		{
      return _val != read_index<B>::read(b);
		}

    template<typename B>
		bool operator <=(const B& b) const
		{
      return _val <= read_index<B>::read(b);
		}

    template<typename B>
		bool operator >=(const B& b) const
		{
      return _val >= read_index<B>::read(b);
		}

    template<typename B>
		bool operator <(const B& b) const
		{
      return _val < read_index<B>::read(b);
		}

    template<typename B>
		bool operator >(const B& b) const
		{
      return _val > read_index<B>::read(b);
		}
  };

  /** @} */

  namespace core
  {
    /** @cond PRIVATE */
    template<typename ... A>
    struct static_col_extents;
    /** @endcond */

    /**
     * Reads the static extents from random event types (i.e. random variables)
     * and returns the size and extents for a given list of these types
     * thus allowing to determine the row and columns extents of the individual
     * variables
     */
    template<typename B, typename ... C>
    struct static_col_extents<B, given, C...>
    {
      typedef std::tuple<int> extents_type;
      static extents_type extents()
      {
        return std::make_tuple(B::extent());
      }

      static int size()
      {
        return B::extent();
      }
    };

    template<typename B>
    struct static_col_extents<B>
    {
      typedef std::tuple<int> extents_type;
      static extents_type extents()
      {
        return std::make_tuple(B::extent());
      }

      static int size()
      {
        return B::extent();
      }
    };

    template<typename A, typename ... C>
    struct static_col_extents<A, C...>
    {
      typedef decltype(
          util::tuple::cons(A::extent(),
              static_col_extents<C...>::extents())) extents_type;

      static extents_type extents()
      {
        return util::tuple::cons(A::extent(),
            static_col_extents<C...>::extents());
      }

      static int size()
      {
        return util::tuple::fold([] (int a, int b)
        -> int
        {
          return a*b;
        }, 1, extents());
      }
    };

    /** @cond PRIVATE */
    template<typename ... A>
    struct static_row_extents;
    /** @endcond */

    template<>
    struct static_row_extents<>
    {
      typedef std::tuple<> extents_type;
      static extents_type extents()
      {
        return std::make_tuple();
      }

      static int size()
      {
        return 1;
      }
    };

    /**
     * Reads the static extents from random event types (i.e. random variables)
     * and returns the size and extents for a given list of these types
     * thus allowing to determine the row and columns extents of the individual
     * variables
     */
    template<typename B>
    struct static_row_extents<given, B>
    {
      typedef std::tuple<int> extents_type;
      static extents_type extents()
      {
        return std::make_tuple(B::extent());
      }

      static int size()
      {
        return B::extent();
      }
    };

    template<typename B, typename ... C>
    struct static_row_extents<given, B, C...>
    {
      typedef decltype(
          util::tuple::cons(B::extent(),
              static_row_extents<given,C...>::extents())) extents_type;

      static extents_type extents()
      {
        return util::tuple::cons(B::extent(),
            static_row_extents<given, C...>::extents());
      }

      static int size()
      {
        return util::tuple::fold([] (int a, int b)
        -> int
        {
          return a*b;
        }, 1, extents());
      }
    };

    template<typename A, typename ... C>
    struct static_row_extents<A, C...>
    {
      typedef typename static_row_extents<C...>::extents_type extents_type;
      static extents_type extents()
      {
        return static_row_extents<C...>::extents();
      }

      static int size()
      {
        return static_row_extents<C...>::size();
      }
    };

    /** @cond PRIVATE */
    template<typename ... T>
    struct index_iterator;
    /** @endcond */

    /**
     * @brief Used to iterate a function over all indices given specific extents
     */
    template<typename Head, typename ...Tail, template<typename ...> class T>
    struct index_iterator<T<Head, Tail...>>
    {
      template<typename F, typename... Index,
      typename ExtentHead, typename... ExtentTail,
      template<typename ...> class E,
      template <typename...> class I = std::tuple>
      static void apply_all(F&& f,
          const I<Index...>& idx,
          const E<ExtentHead, ExtentTail...>& extents)
      {
        unsigned headExtent = read_index<ExtentHead>::read(std::get<0>(extents));
        for(unsigned i=0; i<headExtent; ++i)
        {
          // Append the current index to the index tuple
        	I<Index..., Head> new_idx = util::tuple::append(idx, Head(i));
        	// And remove the current extent from the extent tuple
          E<ExtentTail...> extent_tail = util::tuple::tail(extents);
          // Recursively call apply
          index_iterator<T<Tail...>>::apply_all(
              f,
              new_idx,
              extent_tail
          );
        }
      }

    };

    template<template<typename ...> class T>
    struct index_iterator<T<>>
    {
      /** @cond PRIVATE */
      template<typename S>
      struct apply_impl;
      /** @endcond */

      template<size_t ...IndexIterator>
      struct apply_impl<util::compile_time_list::integer_list<IndexIterator...>>
      {
        template<typename ...Index, typename F,
        template <typename...> class I = std::tuple>
        static void apply(F && f, const I<Index...>& idx_tpl)
        {
          f(std::forward<const Index>(std::get<IndexIterator>(idx_tpl))...);
        }
      };

      template<typename F, typename... Index,
      typename E,
      template <typename...> class I = std::tuple>
      static void apply_all(F && f, const I<Index...>& idx, const E& extents)
      {
        apply_impl<typename util::compile_time_list::iota_0<sizeof...(Index)>::type>::apply(f, idx);
      }
    };

    /** @cond PRIVATE */
    template<typename ... T>
    struct index_iterator_reverse;
    /** @endcond */

    /**
     * @brief Used to iterate a function over all indices given specific extents
     *
     * The inner most loop of index_iterator is the outer most loop here
     */
    template<typename Head, typename ...Tail, template<typename ...> class T>
    struct index_iterator_reverse<T<Head, Tail...>>
    {
      template<typename F, typename... Index,
      typename... Extent,
      template<typename ...> class E,
      template <typename...> class I = std::tuple>
      static void apply_all(F&& f,
          const I<Index...>& idx,
          const E<Extent...>& extents)
      {
        typedef typename std::decay<decltype(util::tuple::last(extents))>::type Last;
        unsigned lastExtent = read_index<Last>::read(util::tuple::last(extents));
        for(unsigned i=0; i<lastExtent; ++i)
        {
        	// Insert the current index at the front
          auto new_idx = util::tuple::cons(i, idx);
          // Remove the current extents from the end of the extents tuple
          auto extent_start = util::tuple::start(extents);
          // Recursively apply
          index_iterator_reverse<T<Tail...>>::apply_all(
              f,
              new_idx,
              extent_start
          );
        }
      }

    };

    template<template<typename ...> class T>
    struct index_iterator_reverse<T<>>
    {
    	/** @cond PRIVATE */
    	template<typename S>
    	struct apply_impl;
    	/** @endcond */

    	template<size_t ...IndexIterator>
    	struct apply_impl<util::compile_time_list::integer_list<IndexIterator...>>
    	{
    		template<typename ...Index, typename F,
    		template <typename...> class I = std::tuple>
    		static void apply(F && f, const I<Index...>& idx_tpl)
    		{
    			f(std::forward<const Index>(std::get<IndexIterator>(idx_tpl))...);
    		}
    	};

    	template<typename F, typename... Index,
    	typename E,
    	template <typename...> class I = std::tuple>
    	static void apply_all(F && f, const I<Index...>& idx, const E& extents)
    	{
    		apply_impl<typename util::compile_time_list::iota_0<sizeof...(Index)>::type>::apply(f, idx);
    	}
    };
  }
}

#endif /* _RANDOMVARIABLE_H_ */
