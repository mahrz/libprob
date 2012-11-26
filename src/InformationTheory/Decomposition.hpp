#ifndef _DECOMPOSITION_H_
#define _DECOMPOSITION_H_

#define PROB_NL_PRECISION 1e-5

#include <nlopt.hpp>

/**
 * @file Decomposition	.hpp
 *
 * @brief Different decomposition/redundancy measures and helper functions
 * @todo This is experimental and not tested
 *
 */

namespace prob
{

  namespace it
  {
    /**
     * @defgroup DECOMP Decomposition
     * @ingroup IT
     *
     * @brief Information decomposition functions (EXPERIMENTAL POSSIBLY BROKEN!!!)
     *
     * @{
     */

    namespace decomp
    {
      namespace core
      {
        double row_sum_one(const std::vector<double> &x,
            std::vector<double> &grad, void *data)
        {
          std::tuple<unsigned, unsigned> * data_int =
              reinterpret_cast<std::tuple<unsigned, unsigned>*>(data);
          unsigned row = std::get<0> (*data_int);
          unsigned cols = std::get<1> (*data_int);

          if (!grad.empty())
          {
            for (unsigned i = row * cols; i < (row + 1) * cols; i++)
              grad[i] = 1;
          }

          double sum = 0;
          for (unsigned i = row * cols; i < (row + 1) * cols; i++)
            sum += x[i];

          return sum - 1;
        }

        double sum_one(const std::vector<double> &x, std::vector<double> &grad,
            void *data)
        {

          if (!grad.empty())
          {
            for (unsigned i = 0; i < x.size(); i++)
              grad[i] = 1;
          }
          double sum = 0;
          for (unsigned i = 0; i < x.size(); i++)
            sum += x[i];

          return sum - 1;
        }

        template<typename DistZ>
        double red_objective(const std::vector<double> &x,
            std::vector<double> &grad, void *data)
        {
          unsigned dim = x.size();
          std::tuple<DistZ const * const, std::vector<DistZ> const * const >* data_tuple =
              reinterpret_cast<std::tuple<DistZ const * const,
                  std::vector<DistZ> const * const >*>(data);

          DistZ q(std::get<1> (*data_tuple)->front());
          q.setZero();

          double sum = 1;
          for (unsigned i = 0; i < dim; i++)
          {
            q += x[i] * (*std::get<1> (*data_tuple))[i];
            sum -= x[i];
          }
          q += sum * (*std::get<1> (*data_tuple))[dim];

          return it::kl_divergence((*std::get<0> (*data_tuple)), q);

        }

        template<typename DistXYgZ, typename DistZ, typename DistZgZ,
            typename DistXgZ, typename DistYgZ, size_t sX, size_t sY, size_t sZ>
        double icm_objective(const std::vector<double> &x,
            std::vector<double> &grad, void *data)
        {
          typedef std::tuple<DistXYgZ const * const, DistZ const * const > data_type;

          data_type * data_tuple = reinterpret_cast<data_type*>(data);

          DistXYgZ const * const dXYgZ = std::get<0> (*data_tuple);
          DistZ const * const dZ = std::get<1> (*data_tuple);

          DistZgZ dZngZ = square(*dZ);

          unsigned dimZ = dZ->cols();

          unsigned i = 0;

          for (unsigned r = 0; r < dimZ; r++)
            for (unsigned c = 0; c < dimZ; c++)
              dZngZ.coeffRef(r, c) = x[i++];

          std::cout << dZngZ << std::endl;

          auto dXYZngZ = join_conditionals(*dXYgZ, dZngZ);

          typedef typename util::compile_time_list::iota_0<sX + sY + sZ>::type index_typeXYZ;
          auto dXYZn = marginalize(uncondition(dXYZngZ, *dZ), index_typeXYZ());
          typedef typename util::compile_time_list::iota_n<sX + sY, sX + sY + sZ>::type index_typeZ;
          DistZ dZn = marginalize(dXYZn, index_typeZ());

          typedef typename util::compile_time_list::iota_n<sX, sX + sY + sZ>::type index_typeYZ;
          auto dYZn = marginalize(dXYZn, index_typeYZ());

          typedef typename util::compile_time_list::join_lists<
              typename util::compile_time_list::iota_0<sX>::type, index_typeZ>::type index_typeXZ;
          auto dXZn = marginalize(dXYZn, index_typeXZ());

          DistXYgZ dXYgZn;
          condition(dXYZn, dZn, dXYgZn);

          DistXgZ dXgZn;
          condition(dXZn, dZn, dXgZn);

          DistYgZ dYgZn;
          condition(dYZn, dZn, dYgZn);

          std::cout << dXYgZn << std::endl;

          double cmi = conditional_mutual_information(dXYgZn, dXgZn, dYgZn,
              dZn);

          std::cout << cmi << std::endl;

          return cmi;
        }

        /** @cond PRIVATE */
        template<typename ...T>
        struct redundancy_impl;

        template<typename ...T>
        struct icm_information_impl;

        template<typename ...T>
        struct minimal_information_impl;
        /** @endcond */

        template<template<typename ...> class V,
        typename ... X,
        typename ... Y,
        typename ... Z,
        typename Scalar>
        struct redundancy_impl<Scalar, V<X...>, V<Y...>, V<Z...>>
        {
          typedef distribution<Scalar, Z..., given, X...> DistZgX;
          typedef distribution<Scalar, Z..., given, Y...> DistZgY;
          typedef distribution<Scalar, Z..., Y...> DistZY;

          typedef distribution<Scalar, X...> DistX;
          typedef distribution<Scalar, Y...> DistY;
          typedef distribution<Scalar, Z...> DistZ;

          typedef Eigen::Matrix<double, Eigen::Dynamic, 1> row_vector;

          static std::vector<DistZ> project(
              const std::vector<DistZ>& source,
              const std::vector<DistZ>& target
          )
          {
            std::vector<DistZ> projected;

            unsigned dim = target.size() - 1;
            std::vector<double> lb(dim,0);
            std::vector<double> ub(dim,1);

            nlopt::opt opt(nlopt::LN_COBYLA, dim);

            opt.set_lower_bounds(lb);
            opt.set_upper_bounds(ub);

            opt.add_inequality_constraint(sum_one, nullptr, PROB_NL_PRECISION);
            opt.set_xtol_rel(PROB_NL_PRECISION);

            for(const DistZ& p : source)
            {
              std::tuple<DistZ const * const, std::vector<DistZ> const * const>
              data = std::make_tuple(&p, &target);

              opt.set_min_objective(red_objective<DistZ>, (void*)&data);

              std::vector<double> x(dim, 0.5);
              double minf;
              opt.optimize(x, minf);

              DistZ q(p);
              q.setZero();

              double sum=1;
              for(unsigned i=0;i<dim;i++)
              {
                q += x[i] * target[i];
                sum -= x[i];
              }
              q += sum * target[dim];

              projected.push_back(q);
            }

            return projected;
          }

          static Scalar redundancy(
              const DistZgX& dZgX,
              const DistZgY& dZgY,
              const DistX& dX,
              const DistY& dY,
              const DistZ& dZ)
          {
            std::vector<DistZ> vdZgX;
            std::vector<DistZ> vdZgY;

            std::vector<DistZ> projected_vdZgX;
            std::vector<DistZ> projected_vdZgY;

            // Extract the conditional distributions
            dZgX.each_conditional_index(
                [&vdZgX, &dZgX] (const X&... x)
                {
                  vdZgX.push_back(dZgX.posterior_distribution(x...));
                });

            dZgY.each_conditional_index(
                [&vdZgY, &dZgY] (const Y&... y)
                {
                  vdZgY.push_back(dZgY.posterior_distribution(y...));
                });

            projected_vdZgX = project(vdZgX, vdZgY);
            projected_vdZgY = project(vdZgY, vdZgX);

            Scalar red_x(0), red_y(0);
            typename std::vector<DistZ>::iterator pX = projected_vdZgX.begin();
            typename std::vector<DistZ>::iterator pY = projected_vdZgY.begin();

            dZgX.each_conditional_index(
                [&pX, &dZgX, &dX, &dZ, &red_x] (const X&... x)
                {
                  if(dX(x...) <= 1e-18)
                  {
                    pX++;
                    return;
                  }

                  DistZ dZgx = dZgX.posterior_distribution(x...);
                  Scalar dx = dX(x...);

                  dZ.each_index(
                      [&pX, &dZgx, &dx, &dZ, &red_x] (const Z&... z)
                      {
                        if(dZgx(z...) >= 1e-18)
                        red_x += dZgx(z...) * dx *
                        log2_fraction((*pX)(z...),dZ(z...));
                      }
                  );

                  pX++;
                });

            dZgY.each_conditional_index(
                [&pY, &dZgY, &dY, &dZ, &red_y] (const Y&... y)
                {
                  if(dY(y...) <= 1e-18)
                  {
                    pY++;
                    return;
                  }

                  DistZ dZgy = dZgY.posterior_distribution(y...);
                  Scalar dy = dY(y...);

                  dZ.each_index(
                      [&pY, &dZgy, &dy, &dZ, &red_y] (const Z&... z)
                      {
                        if(dZgy(z...) >= 1e-18)
                        red_y += dZgy(z...) * dy *
                        log2_fraction((*pY)(z...),dZ(z...));
                      }
                  );

                  pY++;
                });

            return std::min(red_y, red_x);
          }
        };

        template<template<typename ...> class V,
        typename ... X,
        typename ... Y,
        typename ... Z,
        typename Scalar>
        struct icm_information_impl<Scalar, V<X...>, V<Y...>, V<Z...>>
        {
          typedef distribution<Scalar, X... , Y..., given, Z...> DistXYgZ;
          typedef distribution<Scalar, Z..., given, Z...> DistZgZ;
          typedef distribution<Scalar, X..., given, Z...> DistXgZ;
          typedef distribution<Scalar, Y..., given, Z...> DistYgZ;

          typedef distribution<Scalar, X...> DistX;
          typedef distribution<Scalar, Y...> DistY;
          typedef distribution<Scalar, Z...> DistZ;

          static Scalar icm_information(
              const DistXYgZ& dXYgZ,
              const DistZ& dZ)
          {

            unsigned dimZ = dZ.cols();
            unsigned dimZgZ = dimZ * dimZ;
            std::vector<std::tuple<unsigned,unsigned>> rc;
            std::vector<double> lb(dimZgZ,0);
            std::vector<double> ub(dimZgZ,1);

            nlopt::opt opt(nlopt::LN_COBYLA, dimZgZ);

            opt.set_lower_bounds(lb);
            opt.set_upper_bounds(ub);

            for(unsigned z = 0; z < dimZ; ++z)
            {
              rc.push_back(std::make_tuple(z,dimZ));
              opt.add_equality_constraint(row_sum_one, &rc[z], PROB_NL_PRECISION);
          }

          opt.set_xtol_rel(PROB_NL_PRECISION);

          std::tuple<DistXYgZ const * const,
          DistZ const * const> data =
          std::make_tuple(&dXYgZ, &dZ);

          opt.set_min_objective(icm_objective<
              DistXYgZ,
              DistZ,
              DistZgZ,
              DistXgZ,
              DistYgZ,
              sizeof...(X),
              sizeof...(Y),
              sizeof...(Z)>, (void*)&data);

          std::vector<double> x(dimZgZ, 1/dimZ);
          double minf;
          opt.optimize(x, minf);

          return Scalar(minf);
        }

      };

        template<template<typename ...
      > class V,
      typename ... X,
      typename ... Y,
      typename ... Z,
      typename Scalar>
      struct minimal_information_impl<Scalar, V<X...>, V<Y...>, V<Z...>>
      {
        typedef distribution<Scalar, X..., given, Z...> DistXgZ;
        typedef distribution<Scalar, Y..., given, Z...> DistYgZ;

        typedef distribution<Scalar, X...> DistX;
        typedef distribution<Scalar, Y...> DistY;
        typedef distribution<Scalar, Z...> DistZ;

        static Scalar minimal_information(
            const DistXgZ& dXgZ,
            const DistYgZ& dYgZ,
            const DistZ& dZ)
        {
          Scalar min_info(0);

          typedef typename util::compile_time_list::iota_0<sizeof...(X)>::type index_typeX;
        DistX dX(marginalize(uncondition(dXgZ,dZ), index_typeX()));

        typedef typename util::compile_time_list::iota_0<sizeof...(Y)>::type index_typeY;
        DistY dY(marginalize(uncondition(dYgZ,dZ), index_typeY()));

        dZ.each_index([&] (const Z&... z)
            {
              Scalar x_info(0);
              Scalar y_info(0);

              DistX dXgz = dXgZ.posterior_distribution(z...);
              DistY dYgz = dYgZ.posterior_distribution(z...);

              dX.each_index([&] (const X&... x)
                  {
                    x_info += xlogy( dXgz(x...), dXgz(x...) / dX(x...));
                  });

              x_info /= log_of_2<Scalar>();

              dY.each_index([&] (const Y&... y)
                  {
                    y_info += xlogy( dYgz(y...), dYgz(y...) / dY(y...));
                  });

              y_info /= log_of_2<Scalar>();

              min_info += dZ(z...) * std::min(x_info, y_info);
            });

        return min_info;
      }

    };

        }

      /**
       * Calculate the redundant information defined as
       * @f[ I_{ \operatorname{red}}(Z;X,Y) = \min \{ I^\pi_Z(X \searrow Y), I^\pi_Z(Y \searrow X) \}  @f]
       * by <a href="http://arxiv.org/abs/1207.2080">Harder, Salge and Polani</a>.
       * This value can be used to decompose the mutual information @f$I(Z;X,Y)@f$ by
       * using the redundant information as the redundancy term in the decomposition
       * of Williams and Beer. See the linked article for a detailed account.
       *
       * @param dZgX @f$p(z|x)@f$
       * @param dZgY @f$p(z|y)@f$
       * @param dX @f$p(x)@f$
       * @param dY @f$p(y)@f$
       * @param dZ @f$p(z)@f$
       * @return @f$ I_{ \operatorname{red}}(Z;X,Y) @f$ as a Scalar
       */
        template<typename DistZgX, typename DistZgY, typename DistX,
        typename DistY, typename DistZ>
        typename DistZgX::scalar redundancy(const DistZgX& dZgX,
            const DistZgY& dZgY, const DistX& dX, const DistY& dY,
            const DistZ& dZ)
        {
          return core::redundancy_impl<typename DistZ::scalar,
          typename DistX::posterior_type, typename DistY::posterior_type,
          typename DistZgY::posterior_type>::redundancy(dZgX, dZgY, dX, dY,
              dZ);
        }

        /** Calculate the intrinsic conditional mutual information defined as
         *  @f[ I_{ \operatorname{icm}}(Z;X,Y) = \min_{p(z'|z)} I(X;Y|Z')  @f]
         *  with @f$|Z'|=|Z|@f$ by
         *  <a href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.26.3569">Maurer and Wolf</a>.
         *  It has also been introduced as measure of unique information for
         *  the decomposition of the mutual information @f$I(Y;Z,X)@f$ and
         *  @f$I(X;Z,Y)@f$ (what @f$ X @f$ uniquely contains about  @f$ Y @f$ and vice versa) by
         *  <a href="http://arxiv.org/abs/1205.4265">Griffith and Koch</a>.
         *
         * @param dXYgZ @f$p(x,y|z)@f$
         * @param dZ @f$p(z)@f$
         * @param dX Dummy argument to determine the type list of X
         * @param dY Dummy argument to determine the type list of Y
         * @return @f$I_{ \operatorname{icm}}(Z;X,Y) @f$ as a Scalar
         */
      template<typename DistXYgZ, typename DistZ, typename DistX, typename DistY>
      typename DistZ::scalar icm_information(const DistXYgZ& dXYgZ,
          const DistZ& dZ, const DistX& dX, const DistY& dY)
      {
        return core::icm_information_impl<typename DistZ::scalar,
            typename DistX::posterior_type, typename DistY::posterior_type,
            typename DistZ::posterior_type>::icm_information(dXYgZ, dZ);
      }

      /**
       * Calculate the minimal information defined as
       * @f[ I_\min(Z;\{X,Y\}) = \sum_z p(z) \min_{A\in\{X,Y\}} I(Z=z;A)  @f]
       * by <a href="http://arxiv.org/abs/1004.2515">Williams and Beer</a>.
       * This can then be used to decompose the mutual information @f$I(Z;X,Y)@f$ by
       * using the minimal information as the redundancy term. See the linked article
       * for a detailed account.
       * @param dXgZ @f$p(x|z)@f$
       * @param dYgZ @f$p(y|z)@f$
       * @param dZ @f$p(z)@f$
       * @return @f$I_\min(Z;\{X,Y\})@f$ as a Scalar
       */
      template<typename DistXgZ, typename DistYgZ, typename DistZ>
      typename DistZ::scalar minimal_information(const DistXgZ& dXgZ,
          const DistYgZ& dYgZ, const DistZ& dZ)
      {
        return core::minimal_information_impl<typename DistZ::scalar,
            typename DistXgZ::posterior_type, typename DistYgZ::posterior_type,
            typename DistZ::posterior_type>::minimal_information(dXgZ, dYgZ, dZ);
      }

    }

    /** @} */
  }
}
#endif /* _DECOMPOSITION_H_ */
