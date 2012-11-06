#ifndef _INITIALIZERS_H_
#define _INITIALIZERS_H_

#include <random>

/**
 * @file Initializers.hpp
 *
 * @brief Initialize distributions
 *
 */

namespace prob
{
  /**
   * @defgroup INIT Initializers
   * @ingroup DIST
   *
   * @brief Some initializers for discrete probability distributions
   *
   * @{
   */

  /** Initialize distributions */
  namespace init
  {
    /**
     * Initializes a \f$n\f$-dimensional probability distribution
     * with the probability \f$ \frac{1}{n}\f$ in each event.
     * @param d Probability Distribution
     */
    template<typename Scalar, typename ...T>
    void uniform(distribution<Scalar, T...>& d)
    {
      int cols = d.cols();
      d.setConstant(Scalar(1.0 / cols));
    }

    /**
     * Initializes a probability distribution randomly with
     * random numbers supplied by the random number generator
     * (e.g. std::mt19937).
     * @param d
     * @param rng
     */
    template<typename Scalar, typename R, typename ... T>
    void random(distribution<Scalar, T...>& d, R& rng)
    {
      std::uniform_real_distribution<> unit_interval(0, 1);
      for (int i = 0; i < d.rows(); ++i)
        for (int j = 0; j < d.cols(); ++j)
          d.coeffRef(i, j) = unit_interval(rng);
      d.normalize();
    }
  }

  /** @} */
}

#endif /* _INITIALIZERS_H_ */
