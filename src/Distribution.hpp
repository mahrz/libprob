#ifndef _DISTRIBUTION_H_
#define _DISTRIBUTION_H_

#include "prob.hpp"
#include <utility>

/**
 * @file Distribution.hpp
 *
 * @brief The main class of the library providing discrete probability distributions
 *
 */

/**
 * For Version v0.1
 *
 * @todo Write unit tests
 * @todo Main page and examples
 * @todo License & Readme
 * @todo Externalize repository
 * @todo Github pages
 *
 * Important but still not in v0.1
 * @todo More unit tests for the internals
 *
 * New Features
 * @todo Information Decomposition (Reconstruction functions)
 * @todo Advanced Information Theoretic Functions
 * @todo Alternative Backends
 * @todo Contexts
 * @todo Advanced Initializers
 *
 * Should be done but not neccessary
 * @todo Staticly assert equality of types for it function (gives ct error anyway)
 * @todo Use map/map conditionals in IT
 * @todo The whole rvalue reference & template thing is still a bit of a mess (check that everything works, and is efficient)
 * @todo Check for static/dynamic consistency (should give errors anyway, but make it clearer)
 *
 * Would be nice to have
 * @todo Better summary of distributions (especially for conditional distributions)
 */

namespace prob
{
	namespace core
	{
		/** @cond PRIVATE */
		template<typename T>
		struct type_printer;

		template<typename ...T>
		struct index_printer;

		template<bool>
		struct conditional_case;

		/** @endcond */

		/**
		 * Used to print the names of random variable types
		 * declared with the #RVAR(x) macro (see \ref RVAR).
		 */
		template<template<typename ...> class I>
		struct type_printer<I<>>
		{
			static void print_types(std::ostream& out)
			{
			}
		};

		/**
		 * Used to print the names of random variable types
		 * declared with the #RVAR(x) macro (see \ref RVAR).
		 */
		template<template<typename ...> class I, typename Head>
		struct type_printer<I<Head>>
		{
			static void print_types(std::ostream& out)
			{
				out << Head::label();
			}
		};

		/**
		 * Used to print the names of random variable types
		 * declared with the #RVAR(x) macro (see \ref RVAR).
		 */
		template<template<typename ...> class I, typename Head, typename ...T>
		struct type_printer<I<Head, T...>>
		{
			static void print_types(std::ostream& out)
			{
				out << Head::label() << " ";
				type_printer<I<T...>>::print_types(out);
			}
		};


		/**
		 * Used to print the indices represented by a set of random events
		 * (which are instances of random variable types).
		 */
		template<typename Head, typename ...T>
		struct index_printer<Head, T...>
		{
			static void print_index(std::ostream& out, const Head& h, const T&... t)
			{
				out << prob::read_index<Head>::read(h) << " ";
				index_printer<T...>::print_index(out, t...);
			}
		};

		/**
		 * Used to print the indices represented by a set of random events
		 * (which are instances of random variable types).
		 */
		template<typename ...T>
		struct index_printer<prob::given, T...>
		{
			static void print_index(std::ostream& out, const prob::given& g,
					const T&... t)
			{
				out << "| ";
				index_printer<T...>::print_index(out, t...);
			}
		};

		/**
		 * Used to print the indices represented by a set of random events
		 * (which are instances of random variable types).
		 */
		template<typename Head>
		struct index_printer<Head>
		{
			static void print_index(std::ostream& out, const Head& h)
			{
				out << prob::read_index<Head>::read(h);
			}
		};

		/** @cond PRIVATE */
		template<typename ...T>
		struct tuple_read_impl;

		template<typename H, typename ...T>
		struct tuple_read_impl<H, T...>
		{
			static std::tuple<H,T...>  read(std::istream& in)
			{
				int h;
				in >> h;
				return cons(H(h), tuple_read_impl<T...>::read(in));
			}
		};

		template<typename ...T>
		struct tuple_read_impl<given, T...>
		{
			static std::tuple<given,T...>  read(std::istream& in)
			{
				in.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
				return cons(given(1), tuple_read_impl<T...>::read(in));
			}
		};

		template<typename H>
		struct tuple_read_impl<H>
		{
			static std::tuple<H>  read(std::istream& in)
			{
				int h;
				in >> h;
				return std::make_tuple<H>(H(h));
			}
		};

		/** @endcond */

		/** @brief Read a space seperated integer tuple from the input stream */
		template<typename ... T>
		std::tuple<T...> tuple_read(std::istream& in)
		{
			return tuple_read_impl<T...>::read(in);
		}

		/**
		 * Encapsulates different functionality between conditional and non-conditional
		 * distributions.
		 */
		template<>
		struct conditional_case<true>
		{
			template<typename Origin, typename F>
			void each_index(const Origin& o, F f) const
			{
				auto extents = util::tuple::concat(
						util::tuple::append(o._col_extents, 1), o._row_extents);
				core::index_iterator<typename Origin::expanded_type>::apply_all(f,
						std::make_tuple(), extents);
			}

			template<typename Origin, typename F>
			void each_index_reverse(const Origin& o, F f) const
			{
				auto extents = util::tuple::concat(
						util::tuple::append(o._col_extents, 1), o._row_extents);
				core::index_iterator_reverse<typename Origin::expanded_type>::apply_all(
						f, std::make_tuple(), extents);
			}

			template<typename Origin, typename F>
			void each_conditional_index(const Origin& o, F f) const
			{
				core::index_iterator<typename Origin::conditional_type>::apply_all(f,
						std::make_tuple(), o._row_extents);
			}
		};

		/**
		 * @brief Encapsulates different functionality between conditional and non-conditional
		 * distributions.
		 */
		template<>
		struct conditional_case<false>
		{
			template<typename Origin, typename F>
			void each_index(const Origin& o, F f) const
			{
				core::index_iterator<typename Origin::posterior_type>::apply_all(f,
						std::make_tuple(), o._col_extents);
			}

			template<typename Origin, typename F>
			void each_index_reverse(const Origin& o, F f) const
			{
				core::index_iterator_reverse<typename Origin::posterior_type>::apply_all(
						f, std::make_tuple(), o._col_extents);
			}

			template<typename Origin, typename F>
			void each_conditional_index(const Origin& o, F f) const
			{
				static_assert(true, "Not a conditional distribution");
			}
		};

		auto element_index_accu(const std::tuple<random_event, random_event>& t,
				std::tuple<int, int> u)
		-> std::tuple<int,int>
		{
			int lastExtent = std::get<1> (u);
			int curExtent = read_index<random_event>::read(std::get<1> (t));
			int lastIndex = std::get<0> (u);
			int curIndex = read_index<random_event>::read(std::get<0> (t));

			// Assertion on any out of bounds access
			assert(curIndex >= 0 && curIndex < curExtent);

			return std::make_tuple(lastExtent * curIndex + lastIndex,
					lastExtent * curExtent);
		}
	}

	/**
	 * @defgroup DIST Probability Distributions
	 *
	 * @brief Discrete probability distributions
	 *
	 * @{
	 */

	/**
	 * @brief A discrete probability distribution backed by a dense
	 * <a href="http://eigen.tuxfamily.org/">Eigen</a> matrix.
	 *
	 * For the purpose of the documentation we assume
	 * @code T... = X..., given, Y... @endcode (for a conditional
	 * distribution) or @code T... = X...  @endcode
	 * for a non conditional distribution, where all X and Y are random variable
	 * types (see \ref RVAR).
	 *
	 * The distribution is here denoted by @f$ p(X... | Y...) @f$ and its probabilities
	 * by @f$ p(x... | y...) @f$. The posterior distributions now are @f$ p(X... | y...) @f$
	 * also denoted by @f$ p(\cdot | y...) @f$.
	 *
	 * @tparam Scalar The scalar type of the probabilities. All Scalars that
	 * work for Eigen matrices are supported, some functions however assume
	 * floating point Scalars or specifically double or to/from double convertible.
	 * Scalars.
	 * @tparam T... The type list of the random variable types.
	 */
	template<typename Scalar, typename ... T>
	class distribution: public Eigen::Matrix<Scalar,
	core::splitter<T...>::conditional_type::eigen_size,
	core::splitter<T...>::posterior_type::eigen_size>
	{
	public:
		/** @brief Access to the scalar type. */
		typedef Scalar scalar;

		/**
		 * @brief Eigen base matrix type
		 *
		 * The base type of the class, a dynamically sized Eigen::Matrix (matrix or vector). */
		typedef typename Eigen::Matrix<Scalar,
				core::splitter<T...>::conditional_type::eigen_size,
				core::splitter<T...>::posterior_type::eigen_size> matrix_type;

		/**
		 * @brief Posterior distribution matrix base type
		 *
		 * The base type of the posterior distributions i.e @f$ p(\cdot | y...) @f$. */
		typedef typename Eigen::Matrix<Scalar, 1,
				core::splitter<T...>::posterior_type::eigen_size> posterior_matrix_type;

		/**
		 * The base type of the conditional matrix. This type does not represent a
		 * concept from probability theory and is mainly used internally.
		 */
		typedef typename Eigen::Matrix<Scalar,
				core::splitter<T...>::posterior_type::eigen_size, 1> conditional_matrix_type;

		/**
		 * @brief Row extents/indices tuple type
		 *
		 * The tuple type (std::tuple) of the row extents (when denoting the size
		 * of the conditional variables) or row indices when denoting a specific event
		 * (of the conditional variables).
		 */
		typedef typename core::splitter<T...>::conditional_type::index_type row_type;

		/**
		 * @brief Column extents/indices tuple type
		 *
		 * The tuple type (std::tuple) of the column extents (when denoting the size
		 * of the posterior variables) or columns indices when denoting a specific event
		 * (of the posterior variables).
		 */
		typedef typename core::splitter<T...>::posterior_type::index_type col_type;


		/**
		 * @brief Conditional variables type type
		 *
		 * The \ref core::vars type of the conditional variables. Gives access to the splitted
		 * conditionals.
		 */
		typedef typename core::splitter<T...>::conditional_type conditional_type;

		/**
		 * @brief Posterior variables type
		 *
		 * The \ref core::vars type of the posterior variables. Gives access to the splitted
		 * posteriors.
		 */
		typedef typename core::splitter<T...>::posterior_type posterior_type;

		/**
		 * @brief Random variables type
		 *
		 * The \ref core::vars type of all random variables. Gives access to the expanded variables
		 * specificly all \ref _given types are expanded.
		 */
		typedef typename core::splitter<T...>::expanded_type expanded_type;

		/** @cond PRIVATE */
		template<typename _T>
		struct type_to_distribution;
		/** @endcond */

		/**
		 * @brief Variable type to distribution type conversion
		 *
		 * Converts a vars<T...> (or any other U<T...>) type to a
		 * distribution<Scalar, T...> type.
		 */
		template<typename ..._T, template<typename ...> class U>
		struct type_to_distribution<U<_T...>>
		{
			typedef distribution<Scalar, _T...> distribution_type;
		};

		/**
		 * @brief Posterior distribution type
		 *
		 * The type of the posterior distributions i.e @f$ p(\cdot | y...) @f$.
		 */
		typedef typename type_to_distribution<posterior_type>::distribution_type
				posterior_distribution_type;
		typedef typename type_to_distribution<conditional_type>::distribution_type
				conditional_distribution_type;

	private:

		static constexpr bool _conditional_distribution =
				core::splitter<T...>::conditional_distribution;

		template <bool> friend struct core::conditional_case;
		core::conditional_case<_conditional_distribution> cased;

		row_type _row_extents;
		col_type _col_extents;

	public:

		/**
		 * @brief Is the distribution a conditional distribution?
		 *
		 * This is a compile time property of the type.
		 */
		static constexpr bool conditional_distribution()
		{
			return _conditional_distribution;
		}

		/**
		 * @brief Load a distribution
		 *
		 * Load a distribution from any std::istream. See here for the
		 * distribution format.
		 */
		static distribution<Scalar, T...> load(std::istream& in)
    {
			// Ignore the headers
			in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

			// Read the extents
			auto full_extents = core::tuple_read<T...>(in);

			auto col_extents = util::tuple::subset(full_extents,
														typename util::compile_time_list::iota_0<
																				core::splitter<T...>::posteriors()>::type());

			auto row_extents = util::tuple::subset(full_extents,
																	typename util::compile_time_list::iota_n<
																							core::splitter<T...>::posteriors()+1,
																							sizeof...(T)>::type());

			distribution<Scalar, T...> dist;
			// In case of a static distribution this should either match or assert
			dist.reshape_dimensions(row_extents, col_extents);
			dist.setZero();

			bool stop = false;

			do
			{
				auto full_index = core::tuple_read<T...>(in);

				auto col_index = util::tuple::subset(full_index,
															typename util::compile_time_list::iota_0<
																					core::splitter<T...>::posteriors()>::type());

				auto row_index = util::tuple::subset(full_index,
																		typename util::compile_time_list::iota_n<
																								core::splitter<T...>::posteriors()+1,
																								sizeof...(T)>::type());

				stop = !(in >> dist.prob_ref_via_tuple(row_index, col_index));
			} while(!stop);

			return dist;
    }

		/**
		 * @brief Default constructor
		 *
		 * For static random variables this is the default constructor that needs to be used.
		 * For dynamicly sized random variables the default constructor creates a matrix with
		 * empty extents and a one-element vector as backed Eigen matrix. Unless reshape_dimensions
		 * is called at least once the matrix should not be used anywhere.
		 *
		 * @todo Some way to assert any illicit usage?
		 */
		distribution() : matrix_type(core::static_row_extents<T...>::size(),
				core::static_col_extents<T...>::size()),
				_row_extents(core::static_row_extents<T...>::extents()),
				_col_extents(core::static_col_extents<T...>::extents())
		{
		}

		/**
		 * @brief Copy constructor
		 */
		distribution(const distribution<Scalar, T...>& other) :
			matrix_type(other),
			_row_extents(other.row_extents()),
			_col_extents(other.col_extents())
		{
		}

		/**
		 * @brief Quasi copy constructor
		 *
		 * Reshape the distribution but copy the backed matrix
		 * (overall dimensions stay fixed).
		 */
		template<typename S>
		distribution(const S& other,
				const row_type& row_extents,
				const col_type& col_extents) :
				matrix_type(other),
				_row_extents(row_extents),
				_col_extents(col_extents)
		{
			unsigned induced_rows =
					util::tuple::fold(
							[] (const random_event &a, int b)
							{ return a._val*b; },
							1,
							row_extents);
			assert(other.rows() == induced_rows);

			unsigned induced_cols =
					util::tuple::fold(
							[] (const random_event &a, int b)
							{ return a._val*b; },
							1,
							col_extents);
			assert(other.cols() == induced_cols);
		}

		/**
		 * @brief Constructor for distributions with dynamically sized random variables
		 *
		 * Use this constructor to for distributions that use dynamically sized
		 * random variables.
		 *
		 * @code
		 * RVAR(X)
		 * RVAR(Y)
		 * ...
		 * distribution<double, X, given, Y> p(X(10) | Y(5));
		 * @endcode
		 *
		 * Leads to a compile time error if any of the random variable types
		 * is declared as static.
		 *
		 * @tparam _T... Implicit template parameter the expanded type of _T... needs
		 * to match the distribution template type T...
		 * @param t... The extents of each random variable
		 */
		template<typename... _T>
		distribution(_T... t) :
		matrix_type(core::splitter<_T...>::rows(std::forward<_T>(t)...),
				core::splitter<_T...>::cols(std::forward<_T>(t)...)),
				_row_extents(core::splitter<_T...>::row_index(std::forward<_T>(t)...)),
				_col_extents(core::splitter<_T...>::col_index(std::forward<_T>(t)...))
		{
			// Type to split the type list of the constructor parameters
			typedef core::splitter<typename std::decay<_T>::type...> local_splitter;

			// Called on a distribution of static random variables? If yes, compile time error
			extents_assert<T...>::assert_static();

			// Also a compile time error if T... (always expanded)
			// does not match _T... (expanded type)
			static_assert(util::traits::are_equivalent<typename local_splitter::expanded_type, T...>::value,
					"Random variable type mismatch");
		}

		distribution& operator=(const distribution &other)
		{
			if(&other == this)
				return *this;

			matrix_type::operator=(other);
			_row_extents = other.row_extents();
			_col_extents = other.col_extents();
			return *this;
		}

		/** @brief Extents of the conditional variables as a tuple */
		row_type conditional_extents() const { return _row_extents; }
		/** @brief Extents of the posterior variables as a tuple */
		col_type posterior_extents() const { return _col_extents; }

		/** @brief Alias for conditional_extents */
		row_type row_extents() const { return _row_extents; }
		/** @brief Alias for posterior_extents */
		col_type col_extents() const { return _col_extents; }


		/** @brief Extent of the i-th conditional variable */
		template<size_t i>
		int conditional_extent()
		{
			return std::get<i>(_row_extents)._val;
		}

		/** @brief Extent of the i-th posterior variable */
		template<size_t i>
		int posterior_extent()
		{
			return std::get<i>(_col_extents)._val;
		}

		/** @brief Alias for conditional_extent */
		template<size_t i>
		int row_extent()
		{
			return std::get<i>(_row_extents)._val;
		}

		/** @brief Alias for posterior_extent */
		template<size_t i>
		int col_extent()
		{
			return std::get<i>(_col_extents)._val;
		}

		/** @brief Reshape the extents of the distribution
		 *
		 * Reshapes a distribution of dynamically sized random variables and produces
		 * an runtime assertion if called with different extents on a distribution with
		 * staticly sized random variables, unless the extents are matching the static extents
		 * in which case this method does nothing.
		 *
		 * <b>WARNING</b> In contrast to reshape, this method does not retain the probability
		 * values.
		 *
		 * @param rows The number of rows of the backing matrix (product of all conditional extents)
		 * @param cols The number of columns of the backing matrix
		 * @param new_row_extents The tuple of the row extents (conditional extents)
		 * @param new_col_extents The tuple of the column extents (posterior extents)
		 */
		void reshape_dimensions(
				const row_type& new_row_extents,
				const col_type& new_col_extents)
		{
			unsigned rows =
					util::tuple::fold(
							[] (const random_event &a, int b)
							{ return a._val*b; },
							1,
							new_row_extents);

			unsigned cols =
					util::tuple::fold(
							[] (const random_event &a, int b)
							{ return a._val*b; },
							1,
							new_col_extents);

			if(extents_assert<T...>::is_static())
			{
				assert(rows == this->rows());
				assert(cols == this->cols());
			}
			else
			{
				this->resize(rows, cols);
				_row_extents = new_row_extents;
				_col_extents = new_col_extents;
			}
		}

		/** @brief Reshape the extents of the distribution
		 *
		 * Reshapes a distribution of dynamically sized random variables, results in a compile time
		 * error if called on staticly sized distributions. If you do not know what your distribution
		 * might be, and there is no actual reshaping involved if static random variables are supplied
		 * use reshape_dimensions.
		 *
		 * The reshaping conserved the probability values if possible. For an extensions, new indices are
		 * set to probability 0 and exisitng indices are retained. The resulting
		 * distribution is then normalized. If all extents are extended then all probabilities are retained
		 * otherwise, values are changed due to renormalization.
		 *
		 * @tparam _T... Implicit template parameter, the expansion of the type list needs to match T...
		 * @param t... New extents
		 */
		template<typename... _T>
		void reshape(_T&&... t)
		{
			// Splits the methods argument type lists
			typedef core::splitter<typename std::decay<_T>::type...> local_splitter;

			// Contrary to reshape_dimensions, this cannot be called on
			// a distribution of static random variable types, i.e. it gives
			// a compile time error.
			extents_assert<T...>::assert_static();

			// Compile time error if T... does not match _T... (expanded type)
			static_assert(util::traits::are_equivalent<typename local_splitter::expanded_type, T...>::value,
					"Random variable type mismatch");

			// Extract the new extents from the arguments
			row_type new_row_extents = core::splitter<_T...>::row_index(std::forward<_T>(t)...);
			col_type new_col_extents = core::splitter<_T...>::col_index(std::forward<_T>(t)...);

			// Calculate the new columns and rows
			int rows = core::splitter<_T...>::rows(std::forward<_T>(t)...);
			int cols = core::splitter<_T...>::cols(std::forward<_T>(t)...);

			// Create a copy of self
			distribution<Scalar, T...> copy(*this);

			// Update the distribution
			_row_extents = new_row_extents;
			_col_extents = new_col_extents;

			this->resize(rows, cols);
			this->setZero();
			// Copy the probabilities from the backup copy with the correct indices
			this->each_index([this, &copy] (T... t) { this->prob_ref(t...) = copy.prob_read_or_zero(t...); });
			// Final normalization
			this->normalize();
		}

		/**
		 * @brief Read a probability value
		 *
		 * Similar to prob_ref unless that in case of an out of bounds access
		 * Scalar(0) is returned.
		 *
		 * @tparam _T... Implicit template parameter, the expansion of the type list needs to match T...
		 * @param t... Index of the probability
		 * @returns Probability or Scalar(0) in case of out of bounds access
		 */
		template<typename... _T>
		Scalar prob_read_or_zero(_T&&... t)
		{
			// Splits the methods argument type lists
			typedef core::splitter<typename std::decay<_T>::type...> local_splitter;

			// Compile time error if T... does not match _T... (expanded type)
			static_assert(util::traits::are_equivalent<
					typename local_splitter::expanded_type, T...>::value,
					"Random variable type mismatch");

			// Get the conditional indices of the arguments
			row_type row_index(
					local_splitter::row_index(
							std::forward<typename std::decay<_T>::type>(t)...));

			// Get the posterior indices of the arguments
			col_type col_index(
					local_splitter::col_index(
							std::forward<typename std::decay<_T>::type>(t)...));

			int row, col;
			bool out_of_bounds=false;

			// Fold accumulator that creates the row, column index from the
			// conditional/posterior indices and extents this check for out of bounds
			// access
			auto f = [&out_of_bounds] (const std::tuple<random_event, random_event>& t,
					std::tuple<int,int> u)
					-> std::tuple<int,int>
			{
				// Get the accumulated index and extent
				int lastExtent = std::get<1>(u);
				int lastIndex = std::get<0>(u);
				// Get the index and extent of the next random variable
				int curExtent = read_index<random_event>::read(std::get<1>(t));
				int curIndex = read_index<random_event>::read(std::get<0>(t));

				// Assertion on any out of bounds access
				if(curIndex < 0 || curIndex >= curExtent)
					out_of_bounds = true;

				// Return the accumulated index and extent
				return std::make_tuple(
						lastExtent * curIndex + lastIndex,
						lastExtent * curExtent);
			};

			// Use f to calculate row and column index
			row = std::get<0>(util::tuple::fold(
					f,
					std::make_tuple(0,1),
					util::tuple::zip(row_index,_row_extents)));

			col = std::get<0>(util::tuple::fold(
					f,
					std::make_tuple(0,1),
					util::tuple::zip(col_index,_col_extents)));

			// Out of bounds reads a 0
			if(out_of_bounds)
				return Scalar(0);

			return matrix_type::operator()(row,col);
		}

		/**
		 * @brief Get a probability reference using two index tuples
		 *
		 * @param row_index The index tuple of the conditional events
		 * @param row_index The index tuple of the posterior events
		 * @returns Reference to the probability
		 */
		Scalar& prob_ref_via_tuple(const row_type& row_index, const col_type& col_index)
		{
			int row, col;

			// Use the default accumulator to calculate row and column index
			row = std::get<0>(util::tuple::fold(
					core::element_index_accu,
					std::make_tuple(0,1),
					util::tuple::zip(row_index,_row_extents)));

			col = std::get<0>(util::tuple::fold(
					core::element_index_accu,
					std::make_tuple(0,1),
					util::tuple::zip(col_index,_col_extents)));

			return matrix_type::operator()(row,col);
		}

		/**
		 * @brief Get a probability using two index tuples
		 *
		 * @param row_index The index tuple of the conditional events
		 * @param row_index The index tuple of the posterior events
		 * @returns Reference to the probability
		 */
		Scalar prob_via_tuple(const row_type& row_index, const col_type& col_index)
		{
			int row, col;

			// Use the default accumulator to calculate row and column index
			row = std::get<0>(util::tuple::fold(
					core::element_index_accu,
					std::make_tuple(0,1),
					util::tuple::zip(row_index,_row_extents)));

			col = std::get<0>(util::tuple::fold(
					core::element_index_accu,
					std::make_tuple(0,1),
					util::tuple::zip(col_index,_col_extents)));

			return matrix_type::operator()(row,col);
		}

		/**
		 * @brief Get a probability reference
		 *
		 * @tparam _T... Implicit template parameter, the expansion of the type list needs to match T...
		 * @param t... Index of the probability
		 * @returns Reference to the probability
		 */
		template<typename... _T>
		Scalar& prob_ref(_T... t)
		{
			// Splits the methods argument type lists
			typedef core::splitter<typename std::decay<_T>::type...> local_splitter;

			// Compile time error if T... does not match _T... (expanded type)
			static_assert(util::traits::are_equivalent<
					typename local_splitter::expanded_type, T...>::value,
					"Random variable type mismatch");

			// Get the conditional indices of the arguments
			row_type row_index(
					local_splitter::row_index(
							std::forward<typename std::decay<_T>::type>(t)...));

			// Get the posterior indices of the arguments
			col_type col_index(
					local_splitter::col_index(
							std::forward<typename std::decay<_T>::type>(t)...));

			int row, col;

			// Use the default accumulator to calculate row and column index
			row = std::get<0>(util::tuple::fold(
					core::element_index_accu,
					std::make_tuple(0,1),
					util::tuple::zip(row_index,_row_extents)));

			col = std::get<0>(util::tuple::fold(
					core::element_index_accu,
					std::make_tuple(0,1),
					util::tuple::zip(col_index,_col_extents)));

			return matrix_type::operator()(row,col);
		}

		/**
		 * @brief Read a probability value
		 *
		 * @tparam _T... Implicit template parameter, the expansion of the type list needs to match T...
		 * @param t... Index of the probability
		 * @returns Probability value
		 */
		template<typename... _T>
		Scalar operator()(_T ... t) const
		{
			// Splits the methods argument type lists
			typedef core::splitter<typename std::decay<_T>::type...> local_splitter;

			// Compile time error if T... does not match _T... (expanded type)
			static_assert(util::traits::are_equivalent<
					typename local_splitter::expanded_type, T...>::value,
					"Random variable type mismatch");

			// Get the conditional indices of the arguments
			row_type row_index(
					local_splitter::row_index(
							std::forward<typename std::decay<_T>::type>(t)...));

			// Get the posterior indices of the arguments
			col_type col_index(
					local_splitter::col_index(
							std::forward<typename std::decay<_T>::type>(t)...));

			int row, col;

			// Use the default accumulator to calculate row and column index
			row = std::get<0>(util::tuple::fold(
					core::element_index_accu,
					std::make_tuple(0,1),
					util::tuple::zip(row_index,_row_extents)));

			col = std::get<0>(util::tuple::fold(
					core::element_index_accu,
					std::make_tuple(0,1),
					util::tuple::zip(col_index,_col_extents)));

			return matrix_type::operator()(row,col);
		}

		/**
		 * @brief Alias for prob_ref
		 */
		template<typename... _T>
		Scalar& operator()(_T&& ... t)
		{
			return this->prob_ref(t...);
		}

		/**
		 * @brief Get a conditioned posterior distribution
		 *
		 * Results in a compiler error if called on a non-conditional distribution.
		 *
		 * @tparam _T... Implicit template parameter, the expansion of the type list needs
		 * to match the posterior variable type list of T... (i.e. _T... = Y... for
		 * T = X..., given, Y...).
		 * @param t... Conditional indices  @f$ y... @f$
		 * @returns Posterior distribution @f$ p(\cdot | y...) @f$
		 */
		template<typename... _T>
		posterior_distribution_type posterior_distribution(_T&& ... t) const
		{
			// Compile time error if T... does not match _T... (expanded type)
			static_assert(util::traits::are_equivalent<conditional_type,
					typename std::decay<_T>::type...>::value,
					"Random variable type mismatch");

			static_assert(_conditional_distribution, "Not a conditional distribution");

			// The row index does not need to be extracted here, we just make a tuple
			row_type row_index = std::make_tuple(t...);

			int row;

			// Again fold index and extents
			row = std::get<0>(util::tuple::fold(
					core::element_index_accu,
					std::make_tuple(0,1),
					util::tuple::zip(row_index,_row_extents)));

			// Return the row of the matrix as the posterior distribution type with
			// corresponding extents
			return posterior_distribution_type(matrix_type::row(row), std::make_tuple<>(), _col_extents);
		}

		/**
		 * @brief Iterate over all variable indices
		 *
		 * This function iterates over each index t... and calls f(t...).
		 * The variable indices are incremented with the first variable being in the
		 * outer most loop.
		 *
		 * Example using a lambda function:
		 *
		 * @code
		 * distribution<double, X, given, Y> p;
		 * p.each_index([&] (const X& x, given, const Y& y)
		 * 							{
		 * 								cout << x << " " << y << " : " << p(x|y);
		 * 							});
		 * @endcode
		 *
		 * @tparam F function type
		 * @param f the function.
		 */
		template<typename F>
		void each_index(F f) const
		{
			// Different implementations depending on whether the distribution is
			// a conditional or not
			cased.each_index(*this,f);
		}

		/**
		 * @brief Iterate over all variable indices with reversed loops
		 *
		 * This function iterates over each index t... and calls f(t...).
		 * The variable indices are incremented with the first variable being in the
		 * inner most loop.
		 *
		 * Example using a lambda function:
		 *
		 * @code
		 * distribution<double, X, given, Y> p;
		 * p.each_index([&] (const X& x, given, const Y& y)
		 * 							{
		 * 								cout << x << " " << y << " : " << p(x|y);
		 * 							});
		 * @endcode
		 *
		 * @tparam F function type
		 * @param f the function.
		 */
		template<typename F>
		void each_index_reverse(F f) const
		{
			// Different implementations depending on whether the distribution is
			// a conditional or not
			cased.each_index_reverse(*this,f);
		}

		/**
		 * @brief Iterate over all conditional variable indices
		 *
		 * This function iterates over each conditional index t... and calls f(t...).
		 * The variable indices are incremented with the first conditional variable being in the
		 * outer most loop.
		 *
		 * @tparam F function type
		 * @param f the function.
		 */
		template<typename F>
		void each_conditional_index(F f) const
		{
			// Different implementations depending on whether the distribution is
			// a conditional or not
			cased.each_conditional_index(*this,f);
		}

		/**
		 * @brief Normalize the distribution
		 */
		void normalize()
		{
			for(unsigned i=0;i<matrix_type::rows();++i)
			{
				Scalar sum_value = matrix_type::row(i).sum();
				if(sum_value > 0)
					matrix_type::row(i).operator/=(sum_value);
			}
		}

		/**
		 * @brief Sum of all probability values in the distribution
		 *
		 * For a normalized distribution this should always return 1 and
		 * n for a normalized conditional distribution where n is the number
		 * of conditional events.
		 */
		Scalar sum()
		{
			return matrix_type::sum();
		};

		/**
		 * @brief Sum of all probability values grouped by conditionals
		 *
		 * The result is not a proper distribution (called on a normalized
		 * conditional distribution the resulting vector should be containing
		 * only ones), however returning it as a distribution is handy for
		 * certain applications.
		 */
		conditional_distribution_type sum_by_conditional()
		{
			auto summed_matrix =  matrix_type::rowwise().sum();
			return conditional_distribution_type(summed_matrix.transpose(), std::make_tuple<>(), _row_extents);
		};

		/**
		 * @brief Apply a function to each probability value
		 *
		 * The function f only gets the probabilities not the indices.
		 * This method returns the mutated distribution.
		 */
		template<typename F>
		distribution& map(F&& f)
		{
			for(unsigned i=0;i<matrix_type::rows();++i)
				for(unsigned j=0;j<matrix_type::cols();++j)
					matrix_type::operator()(i,j) = f(matrix_type::operator()(i,j));

			return *this;
		}

		/**
		 * @brief Apply a function to a copy of each probability value
		 *
		 * The function f only gets the probabilities not the indices. This method
		 * returns a distribution of the same tyoe.
		 */
		template<typename F>
		distribution map_copy(F&& f) const
		{
			matrix_type mapped(matrix_type::rows(),matrix_type::cols());
			for(unsigned i=0;i<matrix_type::rows();++i)
				for(unsigned j=0;j<matrix_type::cols();++j)
					mapped(i,j) = f(matrix_type::operator()(i,j));
			return distribution(mapped, _row_extents, _col_extents);
		}

		/**
		 * @brief Apply a function to each posterior distribution
		 *
		 * The function f only gets the posterior distribution not the indices of
		 * the conditional. This method returns the mutated distribution.
		 */
		template<typename F>
		distribution& map_by_conditional(F&& f)
		{
			for(unsigned i=0;i<matrix_type::rows();++i)
			{
				matrix_type::row(i) =
						f(posterior_distribution_type(matrix_type::row(i),
								std::make_tuple<>(),
								_col_extents));
			}

			return *this;
		}

		/**
		 * @brief Apply a function to a copy of each posterior distribution
		 *
		 * This method returns a distribution of the same type.
		 */
		template<typename F>
		void map_copy_by_conditional(F&& f) const
		{
			matrix_type mapped(matrix_type::rows(), matrix_type::cols());
			for(unsigned i=0;i<matrix_type::rows();++i)
			{
				mapped.row(i) =
						f(posterior_distribution_type(matrix_type::row(i),
								std::make_tuple<>(),
								_col_extents));
			}
			return distribution(mapped, _row_extents, _col_extents);
		}

		/**
		 * @brief A generalized marginalization method
		 *
		 * Let the random variable type list
		 * @f[ T... = X_0,...,X_k, \operatorname{given}, X_{k+1},..., X_{n} @f]
		 * then grouped_map_sum<i_0,...,i_m>, given that @f$ i_j \leq k @f$ and @f$ i_{j+1} \geq k @f$,
		 * applies the function f as follows:
		 * @f[ p(x_{i_0}, ..., x_{i_j}| x_{i_{j+1}}, ..., x_{i_m}) =
		 * \sum_{x_{l_1}, ..., x_{l_{n-m}}, l_t \notin \{i_0, ... i_m\}} f(p(x_* ...))  @f]
		 * where @f$ x_* ... @f$ denotes the merge of the summation indices and the
		 * selection indices @f$ i_0,...,i_m @f$.
		 *
		 */
		template<int... GroupIndices>
		auto grouped_map_sum(std::function<Scalar(Scalar)> f) const  ->
		typename type_to_distribution<
		typename core::indexed_type_selector<expanded_type,
		core::splitter<T...>::posteriors(),
		core::splitter<T...>::conditionals(),
		-1,
		GroupIndices...>::result_type>::distribution_type
		{
			// Are any indices out of range? No permutation of posteriors/conditionals?
			static_assert(core::check_indices<
					core::splitter<T...>::posteriors(),
					core::splitter<T...>::conditionals(),
					0,
					GroupIndices...>::valid(), "Variable index out of range");

			// Allows to select a subset of the types T... given by the
			// group indices
			typedef typename core::indexed_type_selector<
					expanded_type,
					core::splitter<T...>::posteriors(),
					core::splitter<T...>::conditionals(),
					-1,
					GroupIndices...> type_selector;

			// Using the type selector we define the result type
			typedef typename type_to_distribution<typename type_selector::result_type>
			::distribution_type result_type;

			// First we determine the row (conditional) extents
			// of the result type, by selecting a subset
			// of the current row extents

			// The index splitter splits the group indices
			// between posteriors and conditionals
			typename result_type::row_type grouped_row_extents =
					util::tuple::subset(_row_extents,
							typename
							core::index_splitter<row_type,
							core::splitter<T...>::posteriors(),
							core::splitter<T...>::conditionals(),
							GroupIndices...>::row_index_type());

			/*
			// Then calculate the number of rows of the grouped distribution
			int grouped_rows =
					util::tuple::fold(
							[] (const random_event &a, int b)
							{ return a._val*b; },
							1,
							grouped_row_extents);*/

			// Now we determine the column (posterior) extents
			// of the result type, by selecting a subset
			// of the current column extents
			typename result_type::col_type grouped_col_extents =
					util::tuple::subset(_col_extents,
							typename
							core::index_splitter<col_type,
							core::splitter<T...>::posteriors(),
							core::splitter<T...>::conditionals(),
							GroupIndices...>::col_index_type());
			/*
			// And calculate the number of cols of the grouped distribution
			int grouped_cols =
					util::tuple::fold(
							[] (const random_event &a, int b)
							{ return a._val*b; },
							1,
							grouped_col_extents);*/

			result_type grouped_dist;
			grouped_dist.reshape_dimensions(
					grouped_row_extents, grouped_col_extents);

			grouped_dist.setZero();

			// Finally sum over all indices that are not group indices
			// and apply f to each value
			this->each_index(
					[this, &grouped_dist, &f] (const T&... t)
					{
				core::ref<
				typename type_selector::index_type,
				double,result_type, T...>::get(grouped_dist, t...) +=
						f(this->operator()(std::forward<const T&>(t)...));
					});

			return grouped_dist;
		}

		/** @brief grouped_map_sum with f being the identity */
		template<int... GroupIndices>
		auto grouped_sum() const ->
		typename type_to_distribution<
		typename core::indexed_type_selector<expanded_type,
		core::splitter<T...>::posteriors(),
		core::splitter<T...>::conditionals(),
		-1,
		GroupIndices...>::result_type>::distribution_type
		{
			auto id = [] (Scalar i) { return i; };
			return grouped_map_sum<GroupIndices...>(id);
		}

		/** @brief Alias for grouped_sum */
		template<int... GroupIndices>
		auto marginalize() const ->
		typename type_to_distribution<
		typename core::indexed_type_selector<expanded_type,
		core::splitter<T...>::posteriors(),
		core::splitter<T...>::conditionals(),
		-1,
		GroupIndices...>::result_type>::distribution_type
		{
			auto id = [] (Scalar i) { return i; };
			return grouped_map_sum<GroupIndices...>(id);
		}

		/** @brief Returns a histogram as ASCII art in the given dimensions */
		std::string histogram(unsigned width, unsigned height)
		{
			Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> h(width,height);
			std::stringstream s;

			double max = this->maxCoeff();

			for(unsigned r=0;r<this->rows();r++)
			{
				h.setZero();

				for(unsigned c=0;c<this->cols();c++)
				{
					unsigned hc = (unsigned)((double)(c) / (double) (this->cols()) * height);
					unsigned hr = (unsigned)(this->coeffRef(r,c) / max * (width-1));
					h(hr,hc) = 1;
				}

				for(unsigned hc=0;hc<height;hc++)
				{
					bool first = true;
					for(unsigned hr=0;hr<width;hr++)
					{
						if(h(hr,hc) > 0)
						{
							s << '+';
							first = false;
						}
						else
						{
							if(first)
								s << '-';
							else
								s << ' ';
						}
					}
					s << std::endl;
				}

				if(r != this->rows() - 1)
					for(unsigned x=0;x<width;x++)
						s << '-';
				s << std::endl;
			}

			return s.str();
		}

		/** @brief A summary of the distribution for logging pursposes */
		std::string summary()
		{
			int histogram_width = 20;
			int histogram_height = this->cols() > 20 ? 20 : this->cols();

			std::stringstream s;

			if(prob::distribution<Scalar, T...>::conditional_distribution())
			{
				s << "Conditional Distribution" << std::endl;
				core::type_printer<posterior_type>::print_types(s);
				s << " | ";
				core::type_printer<conditional_type>::print_types(s);
				s << std::endl;
				s << posterior_extents() << " | " << conditional_extents() << std::endl;
			}
			else
			{
				s << "Distribution" << std::endl;
				core::type_printer<posterior_type>::print_types(s);
				s << std::endl;
				s << posterior_extents() << std::endl;
			}

			s << this->histogram(histogram_width, histogram_height);

			return s.str();
		}
	};

	/**
	 * @}
	 */
}

template<typename Scalar, typename ... T>
std::ostream& operator<<(std::ostream& out, const prob::distribution<Scalar, T...>&dist)
{
	typedef prob::distribution<Scalar, T...> dist_type;

	if(prob::distribution<Scalar, T...>::conditional_distribution())
	{
		out << "Conditional Distribution" << std::endl;
		prob::core::type_printer<typename dist_type::posterior_type>::print_types(out);
		out << " | ";
		prob::core::type_printer<typename dist_type::conditional_type>::print_types(out);
		out << std::endl;
		out << dist.posterior_extents() << " | " << dist.conditional_extents() << std::endl;
	}
	else
	{
		out << "Distribution" << std::endl;
		prob::core::type_printer<typename dist_type::posterior_type>::print_types(out);
		out << std::endl;
		out << dist.posterior_extents() << std::endl;
	}

	dist.each_index_reverse(
			[&dist, &out] (const T&... t)
			{
				if(dist(t...) > PROB_EPSILON)
				{
					prob::core::index_printer<T...>::print_index(out, t...);
					out << " : " << dist(t...) << std::endl;
				}
			});

	return out;
}

#endif /* _DISTRIBUTION_H_ */
