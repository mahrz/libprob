#ifndef _INFORMATIONTHEORY_H_
#define _INFORMATIONTHEORY_H_

/**
 * @file InformationTheory.hpp
 *
 * @brief Basic information theoretic functions
 *
 */

namespace prob
{
  /**
   * @defgroup IT Information Theory
   *
   * @brief Information theoretic functions
   *
   * @{
   */



  /**
   * @brief Information theoretic functions
   */
  namespace it
  {

    /** @brief Alias of Eigen::internal::log for convenience */
    template<typename Scalar>
    inline Scalar log(Scalar x)
    {
      return Eigen::internal::log(x);
    }

    /**
     * @brief log(n/d)
     *
     * @f[ \log \frac{n}{d} @f]
     */
    template<typename Scalar>
    inline double log_fraction(Scalar n, Scalar d)
    {
      return log(n / d);
    }

    /** @brief Logarithm (base 10) of 2
     *
     * @f[ \log 2 @f]
     */
    template<typename Scalar>
    inline Scalar log_of_2()
    {
      return Eigen::internal::log(2);
    }

    /** @brief Logarithm (base 10) of 2
     *
     * @f[ \log 2 @f]
     */
    template<>
    inline double log_of_2()
    {
      return 0.693147181;
    }

    /** @brief Logarithm (base 2) of x
     *
     * @f[ \log_2 x @f]
     */
    template<typename Scalar>
    inline Scalar log2(Scalar x)
    {
      return log(x) / log_of_2<Scalar>();
    }

    /**
     * @brief log_2(n/d)
     *
     * @f[ \log_2 \frac{n}{d} @f]
     */
    template<typename Scalar>
    inline Scalar log2_fraction(Scalar n, Scalar d)
    {
      return log_fraction(n, d) / log_of_2<Scalar>();
    }

    /**
     * @brief x*log(y)
     *
     * @f[ \begin{cases} x \log y & \mbox{if } x > 0 \\ 0 & \mbox{otherwise} \end{cases} @f]
     */
    template<typename Scalar>
    inline Scalar xlogy(Scalar x, Scalar y)
    {
      return x > PROB_EPSILON ? x * log(y) : 0;
    }

    /**
     * @brief x*log(x/y)
     *
     * @f[ \begin{cases} x \log \frac{x}{y} & \mbox{if } x > 0 \\ 0 & \mbox{otherwise} \end{cases} @f]
     */
    template<typename Scalar>
    inline Scalar xlogxovery(Scalar x, Scalar y)
    {
      return x > PROB_EPSILON ? x * log(x / y) : 0;
    }

    namespace core
    {

      /** @cond PRIVATE */
      template<typename ...T>
      struct conditional_entropy_impl;

      template<typename ...T>
      struct mutual_information_impl;

      template<typename ...T>
      struct conditional_mutual_information_impl;

      /** @endcond */

      template<template<typename ...> class V,
      typename Scalar,
      typename ...A,
      typename... B>
      struct conditional_entropy_impl<V<A...>, V<B...>, Scalar>
      {
        static Scalar conditional_entropy(
        		const distribution<Scalar,A..., given,B...>& dAgB,
            const distribution<Scalar,B...>& dB)
        {
          Scalar entropy(0);

          dAgB.each_index([&] (const A&... a, given g, const B&... b)
          		{
          			entropy += xlogy<Scalar>(dAgB(a..., g, b...) *
          					dB(b...), dAgB(a..., g, b...));
              });

          return -entropy / log_of_2<Scalar>();
        }
      };

      template<template<typename ...> class V,
      typename Scalar,
      typename ...A, typename... B>
      struct mutual_information_impl<V<A...>, V<B...>, Scalar>
      {
        static Scalar mutual_information(
        		const distribution<Scalar,A..., given, B...>& dAgB,
            const distribution<Scalar,B...>& dB)
        {
          distribution<Scalar,A..., B...> dAB(uncondition(dAgB, dB));
          typedef typename util::compile_time_list::iota_0<sizeof...(A)>::type index_type;
          distribution<Scalar,A...> dA(marginalize(dAB, index_type()));

          return mutual_information(dAgB,
              dA,
              dB);
        }

        static Scalar mutual_information(
        		const distribution<Scalar,A..., given, B...>& dAgB,
            const distribution<Scalar,A...>& dA,
            const distribution<Scalar,B...>& dB)
        {
          Scalar mi(0);

          dAgB.each_index([&] (const A&... a, given g, const B&... b)
              {
          			mi += xlogy<Scalar>(dAgB(a..., g, b...) * dB(b...),
          														dAgB(a..., g, b...) / dA(a...));
              });

          return mi / log_of_2<Scalar>();
        }
      };

      template<template<typename ...> class V,
      typename Scalar,
      typename ...X, typename ...Y, typename ...Z>
      struct conditional_mutual_information_impl<V<X...>, V<Y...>, V<Z...>, Scalar>
      {

        typedef distribution<Scalar, X..., Y..., given, Z...> DistXYgZ;
        typedef distribution<Scalar, X..., given, Z...> DistXgZ;
        typedef distribution<Scalar, Y..., given, Z...> DistYgZ;
        typedef distribution<Scalar, Z...> DistZ;

        static Scalar conditional_mutual_information(
        		const DistXYgZ& dXYgZ,
            const DistXgZ& dXgZ,
            const DistYgZ& dYgZ,
            const DistZ& dZ)

        {
          Scalar mi(0);

          dXYgZ.each_index([&] (const X&... x, const Y&... y, given g, const Z&... z)
              {
								mi += xlogy<Scalar>(dXYgZ(x..., y..., g, z...) * dZ(z...),
										dXYgZ(x..., y..., g, z...) /
										(dXgZ(x..., g, z...) * dYgZ(y..., g, z...)));
              });

          return mi / log_of_2<Scalar>();
        }
      };
    }

    /**
     * @brief Calculate the entropy of a probability distribution
     *
     * @param dist The distribution
     * @return The entropy of dist in bits
     */
    template<typename Scalar, typename ...T>
    Scalar entropy(const distribution<Scalar, T...>& dist)
    {
      static_assert(!distribution<Scalar, T...>::conditional_distribution(),
          "Cannot calculate entropy of a conditional distribution");

      return -dist.map_copy([] (Scalar v)
          { return xlogy<Scalar>(v,v);}).sum() / log_of_2<Scalar>();
    }

    /**
     * @brief Calculate the conditional entropy between to sets of random variables
     *
     * @param dAgB The conditional distribution @f$ p(a...|b...) @f$
     * @param dB The marginal distribution @f$ p(b...) @f$
     * @return @f$ H(A|B) @f$ in bits
     */
    template<typename DistAgB, typename DistB>
    typename DistAgB::scalar conditional_entropy(const DistAgB& dAgB,
        const DistB& dB)
    {
      return core::conditional_entropy_impl<typename DistAgB::posterior_type,
          typename DistAgB::conditional_type, typename DistAgB::scalar>::conditional_entropy(
              dAgB, dB);
    }

    /**
     * @brief Calculate the mutual information between two sets of random variables
     *
     * @param dAgB The conditional distribution @f$ p(a...|b...) @f$
     * @param dB The marginal distribution @f$ p(b...) @f$
     * @return @f$ I(A;B) @f$ in bits
     */
    template<typename DistAgB, typename DistB>
    typename DistAgB::scalar mutual_information(const DistAgB& dAgB,
        const DistB& dB)
    {
      return core::mutual_information_impl<typename DistAgB::posterior_type,
          typename DistAgB::conditional_type, typename DistAgB::scalar>::mutual_information(
              dAgB, dB);
    }

    /**
     * @brief Calculate the mutual information between two sets of random variables
     *
     * In case the other marginal is already calculated you can save some precious
     * processing time by already supplying dA.
     *
     * @param dAgB The conditional distribution @f$ p(a...|b...) @f$
     * @param dA The marginal distribution @f$ p(a...) @f$
     * @param dB The marginal distribution @f$ p(b...) @f$
     * @return @f$ I(A;B) @f$ in bits
     */
    template<typename DistAgB, typename DistA, typename DistB>
    typename DistAgB::scalar mutual_information(const DistAgB& dAgB,
        const DistA& dA, const DistB& dB)
    {
      return core::mutual_information_impl<typename DistAgB::posterior_type,
          typename DistAgB::conditional_type, typename DistAgB::scalar>::mutual_information(
              dAgB, dA, dB);
    }

    /**
     * @brief Calculate the conditional mutual information between three
     * sets of random variables
     *
     * @param dXYgZ The conditional distribution @f$ p(x...,y...|z...) @f$
     * @param dXgZ The marginal distribution @f$ p(x...|z...) @f$
     * @param dYgZ The marginal distribution @f$ p(y...|z...) @f$
     * @param dZ The marginal distribution @f$ p(z...) @f$
     * @return @f$ I(X;Y|Z) @f$ in bits
     */
    template<typename DistXYgZ,
    typename DistXgZ,
    typename DistYgZ,
    typename DistZ>
    typename DistZ::scalar conditional_mutual_information(const DistXYgZ& dXYgZ,
        const DistXgZ& dXgZ, const DistYgZ& dYgZ, const DistZ& dZ)
    {
      return core::conditional_mutual_information_impl<
          typename DistXgZ::posterior_type, typename DistYgZ::posterior_type,
          typename DistZ::posterior_type, typename DistZ::scalar>::conditional_mutual_information(
              dXYgZ, dXgZ, dYgZ, dZ);
    }

    /**
     * @brief Calculate the Kullback-Leibler divergence between two distributions
     * defined on the same random variables.
     *
     * @param dP @f$ p(x...) @f$
     * @param dQ @f$ q(x...) @f$
     * @return @f$ \operatorname{Div}_{KL}(p(\cdot) || q(\cdot)) @f$ in bits
     */
    template<typename Dist>
    typename Dist::scalar kl_divergence(const Dist& dP, const Dist& dQ)
    {
      assert(dP.size() == dQ.size());

      typename Dist::scalar div = 0;
      for (int x = 0; x < dP.cols(); ++x)
        div += xlogxovery(dP.coeffRef(x), dQ.coeffRef(x));

      return div / log_of_2<typename Dist::scalar>();
    }

    /**
     * @brief Calculate the Jensen-Shannon divergence between two distributions
     * defined on the same random variables.
     *
     * @param dP @f$ p(x...) @f$
     * @param dQ @f$ q(x...) @f$
     * @param pi @f$ \pi @f$
     * @return @f$ \operatorname{Div}^\pi_{JS}(p(\cdot) || q(\cdot)) @f$ in bits
     */
    template<typename Dist>
    typename Dist::scalar js_divergence(const Dist& dP, const Dist& dQ,
        typename Dist::scalar pi = 0.5)
    {
    	assert(dP.size() == dQ.size());

      Dist dM(dP);

      for (int x = 0; x < dP.cols(); ++x)
        dM.coeffRef(x) = pi * dP.coeffRef(x) + (1 - pi) * dQ.coeffRef(x);

      return entropy(dM) - pi * entropy(dP) - (1 - pi) * entropy(dQ);
    }

  }

  /** @} */
}

#endif /* _INFORMATIONTHEORY_H_ */
