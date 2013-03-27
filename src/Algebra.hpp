#ifndef _ALGEBRA_H_
#define _ALGEBRA_H_

#include "prob.hpp"

/**
 * @todo Static assertions on type list inequalities
 * @todo Clean up order of template lists
 * @todo Write tests
 * @todo Return matrix by reference set
 * @todo Clean up typenames, and use typedefs using A..., B... for all types
 * @todo Local comments
 */

namespace prob
{
  namespace core
  {
    /** @cond PRIVATE */
    template<typename ...T>
    struct join_conditionals_impl;

    template<typename A, typename B>
    struct join_impl;

    template<typename ...T>
    struct uncondition_impl;

    template<typename ...T>
    struct partial_uncondition_impl;

    template<typename ...T>
    struct condition_impl;

    template<typename ...T>
    struct condition_conditionals_impl;

    template<typename ...T>
    struct bayes_impl;

    template<template<typename ...> class V,
    typename DistA,
    typename DistB,
    typename Scalar,
    typename ...A,
    typename ...B,
    typename ...C>
    struct join_conditionals_impl<V<C...>, V<A...>, V<B...>, Scalar, DistA, DistB>
    {
      typedef prob::distribution<Scalar, A..., B..., given, C...> return_type;

      static return_type join_conditionals(const DistA& distA, const DistB& distB)
      {
        //int rows = distA.rows();
        //int cols = distA.cols() * distB.cols();
        auto col_extents = util::tuple::concat(distA.col_extents(),distB.col_extents());
        return_type result;
        result.reshape_dimensions(distA.row_extents(), col_extents);

        result.each_index(
          [&] (const A&...  a, const B&... b, given g, const C&... c)
          {
            result.prob_ref(a..., b..., g, c...) =
              distA(a...,g,c...) * distB(b...,g,c...);
          });

        return result;
      }
    };

    template<template<typename ...> class V,
    typename Scalar,
    typename ...A,
    typename ...B>
    struct join_impl<V<Scalar, A...>,V<Scalar, B...>>
    {
      typedef V<Scalar, A..., B...> return_type;

      static return_type join(const V<Scalar, A...>& distA, const V<Scalar, B...>& distB)
      {
        auto col_extents = util::tuple::concat(distA.col_extents(),distB.col_extents());
        return_type result;

        result.reshape_dimensions(std::make_tuple<>(), col_extents);

        result.each_index(
            [&] (const A&... a, const B&... b)
            {
              result.prob_ref(a...,b...) =
              distA(a...) * distB(b...);
            });

        return result;
      }
    };

    template<template<typename ...> class V,
    typename DistAgB, typename DistB, typename Scalar,
    typename ...A, typename... B>
    struct uncondition_impl<V<B...>, V<A...>, Scalar, DistAgB, DistB>
    {
      typedef distribution<Scalar, A..., B...> return_type;

      static return_type uncondition(const DistAgB& distAgB, const DistB& distB)
      {
        //int cols = distAgB.rows()*distAgB.cols();

        auto col_extents = util::tuple::concat(distAgB.col_extents(), distAgB.row_extents());
        return_type result;
        result.reshape_dimensions(distB.row_extents(), col_extents);

        result.each_index(
            [&] (const A&... a, const B&... b)
            {
              given g(0);
              result.prob_ref(a..., b...) =
              distAgB(a...,g,b...) * distB(b...);
            });

        return result;
      }
    };

    template<template<typename ...> class V,
    typename DistAgBC, typename DistBgC, typename Scalar,
    typename ...A, typename... B, typename... C>
    struct partial_uncondition_impl< V<A...>, V<C...>, V<B...>, Scalar, DistAgBC, DistBgC>
    {
      typedef distribution<Scalar, A..., B..., given, C...> return_type;

      static return_type partial_uncondition(const DistAgBC& distAgBC, const DistBgC& distBgC)
      {
        //int rows = distBgC.rows();
        //int cols = distAgBC.cols()*distBgC.cols();

        auto col_extents = util::tuple::concat(distAgBC.col_extents(), distBgC.col_extents());
        return_type result;
        result.reshape_dimensions(distBgC.row_extents(), col_extents);

        result.each_index(
            [&] (const A&... a, const B&... b,
                given g, const C&... c)
            {
              result.prob_ref(a..., b..., g, c...) = distAgBC(a..., g, b..., c...) * distBgC(b..., g, c...);
            });

        return result;
      }
    };

    template<template<typename ...> class V,
    typename DistAB, typename DistB, typename Scalar,
    typename ...A, typename... B>
    struct condition_impl<V<A...>, V<B...>, Scalar, DistAB, DistB>
    {
      typedef distribution<Scalar, A..., given, B...> return_type;

      static void condition(const DistAB& distAB, const DistB& distB,
          return_type& distAgB)
      {

        //unsigned rows = distB.cols();
        //unsigned cols = distAB.cols() / rows;
        auto col_extents = util::tuple::subset(
            distAB.col_extents(),
            typename util::compile_time_list::iota_0<sizeof...(A)>::type());

        distAgB.reshape_dimensions(distB.col_extents(), col_extents);

        distAgB.each_index(
            [&] (const A&... a, given g, const B&... b)
            {
              Scalar v = distB(b...);
              if(v==0)
              distAgB.prob_ref(a..., g, b...) = 0;
              else
              distAgB.prob_ref(a..., g, b...) = distAB(a...,b...) / v;
            });

      }
    };

    template<template<typename ...> class V,
    typename DistABgC, typename DistBgC, typename Scalar,
    typename ...A, typename... B, typename... C>
    struct condition_conditionals_impl<V<A...>, V<B...>, V<C...>, Scalar, DistABgC, DistBgC>
    {
      typedef distribution<Scalar, A..., given, B..., C...> return_type;

      static void condition_conditionals(const DistABgC& distABgC,
          const DistBgC& distBgC,
          return_type& distAgBC)
      {

        //unsigned cols = distABgC.cols() / distBgC.cols();
        //unsigned rows = distBgC.cols() * distBgC.rows();
        auto col_extents = util::tuple::subset(
            distABgC.col_extents(),
            typename util::compile_time_list::iota_0<sizeof...(A)>::type());

        auto row_extents = util::tuple::concat(
            distBgC.col_extents(), distBgC.row_extents());

        distAgBC.reshape_dimensions(row_extents, col_extents);

        distAgBC.each_index(
            [&] (const A&... a, given g,
                const B&... b, const C&... c)
            {
              Scalar v = distBgC(b..., g, c...);
              if(v==0)
              distAgBC.prob_ref(a..., g, b..., c...) = 0;
              else
              distAgBC.prob_ref(a..., g, b..., c...) = distABgC(a...,b..., g, c...) / v;
            });
      }
    };

    template<typename Scalar, typename ...A, typename ... B, template<typename ...
  > class V>
  struct bayes_impl<V<A...>, V<B...>, Scalar>
  {
    typedef distribution<Scalar, B..., given, A...> return_type;
    typedef distribution<Scalar, A..., given, B...> conditional_type;
    typedef distribution<Scalar, A...> marginalA_type;
    typedef distribution<Scalar, B...> marginalB_type;

    static return_type bayes(const conditional_type& distAgB,
        const marginalA_type& distA,
        const marginalB_type& distB)
    {
      //int rows = distAgB.cols();
      //int cols = distAgB.rows();

      return_type result;
      result.reshape_dimensions(distAgB.col_extents(), distAgB.row_extents());

        result.each_index(
          [&] (const B&...  b, given g, const A&... a)
          {
            Scalar v = distA(a...);
            if(v==0)
              result.prob_ref(b..., g, a...) = 0;
            else
              result.prob_ref(b..., g, a...) = distAgB(a...,g,b...) * distB(b...) / v;
          });

        return result;
      }
    };
    /** @endcond */
  }

  /**
   * @defgroup DIST Probability Distributions
   *
   * @{
   */

  /**
   * Returns the joint distribution @f$ p(a..., b...) = p(a...)p(b...)@f$.
   *
   * @param dA @f$ p(a...) @f$ of type \ref distribution<Scalar, A...>
   * @param dB @f$ p(b...) @f$ of type \ref distribution<Scalar, B...>
   * @return @f$ p(a..., b...) @f$ of type \ref distribution<Scalar, A..., B...>
   */
  template<typename DistA, typename DistB>
    auto join(const DistA& dA, const DistB& dB) ->
    decltype(core::join_impl<DistA, DistB>::join(dA,dB))
    {
      return core::join_impl<DistA, DistB>::join(dA,dB);
    }

  /**
   * Returns the conditional joint distribution @f$ p(a..., b...|c...) = p(a...|c...)p(b...|c...)@f$.
   *
   * @param dA @f$ p(a...|c...) @f$ of type \ref distribution<Scalar, A..., \ref given, C...>
   * @param dB @f$ p(b...|c...) @f$ of type \ref distribution<Scalar, B..., \ref given, C...>
   * @return @f$ p(a..., b...|c...) @f$ of type \ref distribution<Scalar, A..., B..., \ref given, C...>
   * @todo Proper naming scheme
   */
  template<typename DistA, typename DistB>
  auto join_conditionals(const DistA& dA, const DistB& dB) ->
  decltype(core::join_conditionals_impl<typename DistA::conditional_type,
      typename DistA::posterior_type,
      typename DistB::posterior_type,
      typename DistA::scalar,DistA, DistB>::join_conditionals(dA,dB))
  {
    return core::join_conditionals_impl<typename DistA::conditional_type,
        typename DistA::posterior_type, typename DistB::posterior_type,
        typename DistA::scalar, DistA, DistB>::join_conditionals(dA, dB);
  }

  /**
   * Returns the joint distribution @f$ p(a..., b...) = p(a...|b...)p(b...)@f$.
   *
   * @param dAgB @f$ p(a...|b...) @f$ of type \ref distribution<Scalar, A..., \ref given, B...>
   * @param dB @f$ p(b...) @f$ of type \ref distribution<Scalar, B...>
   * @return @f$ p(a..., b...) @f$ of type \ref distribution<Scalar, A..., B...>
   */
  template<typename DistAgB, typename DistB>
  auto uncondition(const DistAgB& dAgB, const DistB& dB) ->
  decltype(core::uncondition_impl<
      typename DistAgB::conditional_type,
      typename DistAgB::posterior_type,
      typename DistAgB::scalar,
      DistAgB, DistB>::uncondition(dAgB,dB))
  {
    return core::uncondition_impl<typename DistAgB::conditional_type,
        typename DistAgB::posterior_type, typename DistAgB::scalar, DistAgB,
        DistB>::uncondition(dAgB, dB);
  }

  /**
   * Returns the conditional joint distribution @f$ p(a..., b...|c...) = p(a...|b...,c...)p(b...|c...)@f$.
   *
   * @param dAgBC @f$ p(a...|b...,c...) @f$ of type \ref distribution<Scalar, A..., \ref given, B..., C...>
   * @param dBgC @f$ p(b...|c...) @f$ of type \ref distribution<Scalar, B..., \ref given, C...>
   * @return @f$ p(a..., b...|c...) @f$ of type \ref distribution<Scalar, A..., B..., \ref given, C...>
   */
  template<typename DistAgBC, typename DistBgC>
  auto partial_uncondition(const DistAgBC& dAgBC, DistBgC& dBgC) ->
  decltype(core::partial_uncondition_impl<
      typename DistAgBC::posterior_type,
      typename DistBgC::conditional_type,
      typename DistBgC::posterior_type,
      typename DistAgBC::scalar,
      DistAgBC, DistBgC>::partial_uncondition(dAgBC,dBgC))
  {
    return core::partial_uncondition_impl<typename DistAgBC::posterior_type,
        typename DistBgC::conditional_type, typename DistBgC::posterior_type,
        typename DistAgBC::scalar, DistAgBC, DistBgC>::partial_uncondition(
        dAgBC, dBgC);
  }

  /**
   * Returns the conditional distribution @f$ p(a...| b...) = \frac{p(a...,b...)}{p(b...)}@f$.
   *
   * @todo Case of p(b...) = 0;
   *
   * @param dAB @f$ p(a...,b...) @f$ of type \ref distribution<Scalar, A..., B...>
   * @param dB @f$ p(b...) @f$ of type \ref distribution<Scalar, B...>
   * @return @f$ p(a...| b...) @f$ of type \ref distribution<Scalar, A..., \ref given, B...>
   */
  template<typename DistAB, typename DistB, typename DistAgB>
  void condition(const DistAB& dAB, const DistB& dB, DistAgB& dAgB)
  {
    core::condition_impl<typename DistAgB::posterior_type,
        typename DistAgB::conditional_type, typename DistAB::scalar, DistAB,
        DistB>::condition(dAB, dB, dAgB);
  }

  /**
   * Returns the conditional distribution @f$ p(a...| b..., c...) = \frac{p(a...,b...|c...)}{p(b...|c...)}@f$.
   *
   * @todo Case of p(b...|c...) = 0;
   *
   * @param dABgC @f$ p(a...,b...|c...) @f$ of type \ref distribution<Scalar, A..., B..., \ref given, C...>
   * @param dBgC @f$ p(b...) @f$ of type \ref distribution<Scalar, B..., \ref given, C...>
   * @param dAgBC Reference to @f$ p(a...| b..., c...) @f$ of type \ref distribution<Scalar, A..., \ref given, B..., C...>
   */
  template<typename DistABgC, typename DistBgC, typename DistAgBC>
  void condition_conditionals(const DistABgC& dABgC, const DistBgC& dBgC,
      DistAgBC& dAgBC)
  {
    core::condition_conditionals_impl<typename DistAgBC::posterior_type,
        typename DistBgC::posterior_type, typename DistBgC::conditional_type,
        typename DistABgC::scalar, DistABgC, DistBgC>::condition_conditionals(
        dABgC, dBgC, dAgBC);
  }

  /**
   * Returns the conditional distribution @f$ p(b...| a...) = \frac{p(a...|b...)p(b...)}{p(a...)}@f$.
   *
   * @todo Case of p(a...) = 0;
   *
   * @param dAgB @f$ p(a...,b...) @f$ of type \ref distribution<Scalar, A..., \ref given, B...>
   * @param dA @f$ p(a...) @f$ of type \ref distribution<Scalar, A...>
   * @param dB @f$ p(b...) @f$ of type \ref distribution<Scalar, B...>
   * @return @f$ p(b...| a...) @f$ of type \ref distribution<Scalar, B..., \ref given, A...>
   */
  template<typename DistAgB, typename DistA, typename DistB>
  auto bayes(const DistAgB& dAgB, const DistA& dA, const DistB& dB) ->
  decltype(core::bayes_impl<
      typename DistAgB::posterior_type,
      typename DistAgB::conditional_type,
      typename DistAgB::scalar>::bayes(dAgB, dA, dB))
  {
    return core::bayes_impl<typename DistAgB::posterior_type,
        typename DistAgB::conditional_type, typename DistAgB::scalar>::bayes(
        dAgB, dA, dB);
  }

  /**
   * Returns the marginal distribution denoted by the indices in the
   * index list. The index list is a template parameter list of
   * non-negative integers denoting the positions of random variables
   * that form the marginal distribution.
   *
   * Example:
   * @code
   * distribution<double, A, B, given, C, D> p;
   * distribution<double, A, given, C, D> q;
   *
   * q = marginalize(p, util::compile_time_list::integer_list<0,2,3>());
   * @endcode
   * i.e. the 0th index is A, the 2nd C and the 3rd D. The marginalization can be expressed in terms of a sum as
   * @f[ q(a|c,d) = \sum_b p(a,b|c,d). @f]
   *
   * @see distribution.marginalize.
   *
   * @param d The probability distribution
   * @param index An instance of the index list
   * @return The marginalized distribution
   */
  template<typename Dist, template<size_t...> class I, size_t... Indices>
  auto marginalize(const Dist& d, I<Indices...> index) ->
  decltype(d.template marginalize<Indices...>())
  {
    return d.template marginalize<Indices...>();
  }

	/**
	 * @brief Create a self mapping conditional distribution
	 *
	 * Takes a posterior distribution and resizes the current distribution
	 * so that it is in form of a square matrix. This is not necessary for
	 * staticly sized random variables because the extents can be derived from
	 * the type. For dynamically sized variables consider the following example:
	 *
	 * @code
	 * distribution<double, X> p(X(10));
	 * distribution<double, X, given, X> q = distribution<double, X, given, X>::square(p);
	 * @endcode
	 *
	 * is equivalent to
	 *
	 * @code
	 * distribution<double, X> p(X(10));
	 * distribution<double, X, given, X> q(X(10)|X(10));
	 * @endcode
	 *
	 *
	 * However, there are cases where the actual variable list of p is implicit in which
	 * case the square method is a handy shortcut.
	 *
	 * The posterior variables needs to match the conditional variables. Calling this method
	 * for example on the type distribution<double, X, given, Y> will lead to a compile time
	 * error.
	 *
	 * @param other A sized dynamic distribution of the posterior type (Variables: X...)
	 * @returns An sized uninitialized dynamic distribution (Variables: X..., given, X...)
	 */
  template<typename Scalar, typename ...T>
  distribution<Scalar, T..., given, T...> square(const distribution<Scalar, T...>& other)
  {
  	typedef distribution<Scalar, T..., given, T...> dist_type;
  	typename dist_type::col_type col_extents = other.col_extents();
  	typename dist_type::row_type row_extents = other.col_extents();
  	dist_type dist;

		dist.reshape_dimensions(row_extents, col_extents);

		return dist;
  }

  /** @} */

}

#endif /* _ALGEBRA_H_ */
