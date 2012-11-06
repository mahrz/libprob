#ifndef _SPLITTER_H_
#define _SPLITTER_H_

/**
 * @file Splitter.hpp
 *
 * @brief Splitting of packed type lists
 *
 * A lot of internals, splitting of type lists into
 * conditional and posterior types, dealing with given
 * types as well a reference access using indexed subsets
 * of packed arguments.
 */

namespace prob
{
  namespace core
  {
    /** @cond PRIVATE */
    template<typename ... P>
    struct vars;

    template<typename ... A>
    struct splitter;

    template<size_t P, size_t C, size_t ...I>
    struct check_indices;

    template<typename T, bool>
    struct inject_given;

    template<typename T, size_t, bool>
    struct inject_given_index;

    template<typename Scalar, typename D, typename T, typename I>
    struct tuple_ref_getter;

    template<typename T, size_t P, size_t C, int Last, int ... I>
    struct indexed_type_selector;

    template<typename A, typename B, size_t P, bool Given, size_t Head>
    struct cased_append;

    template<typename T, size_t P, size_t C, int ... I>
    struct index_splitter;
    /** @endcond */

    /** @brief Dimensionality type */
    template<>
    struct vars<>
    {
      static const size_t dim = 0;
      static const int eigen_size = 1;

      typedef std::tuple<> index_type;

    };

    /** @brief Dimensionality type */
    template<typename H, typename ... P>
    struct vars<H, P...>
    {
      static const size_t dim = 1 + vars<P...>::dim;
      static const int eigen_size = Eigen::Dynamic;

      typedef typename util::traits::join<std::tuple, std::tuple<H>,
          typename vars<P...>::index_type>::type index_type;
    };

    /**
     * @brief Splits a type list into conditional and posterior variables
     *
     * The splitter is the core struct for getting access to the posterior
     * and conditional variables from a single packed type list (T...). It
     * also provides methods to convert packed arguments into index tuples
     * and to return the number of columns and rows the backend matrix should
     * have.
     */
    template<>
    struct splitter<>
    {
      static constexpr int cols()
      {
        return 1;
      }

      static constexpr int rows()
      {
        return 1;
      }

      static constexpr int posteriors()
      {
        return 0;
      }

      static constexpr int conditionals()
      {
        return 0;
      }

      typedef vars<> posterior_type;
      typedef vars<> conditional_type;

      typedef vars<> expanded_type;

      static constexpr bool conditional_distribution = false;

      static typename posterior_type::index_type col_index()
      {
        return std::make_tuple<>();
      }

      static typename conditional_type::index_type row_index()
      {
        return std::make_tuple<>();
      }
    };

    /**
     * @brief Splits a type list into conditional and posterior variables
     *
     * The splitter is the core struct for getting access to the posterior
     * and conditional variables from a single packed type list (T...). It
     * also provides methods to convert packed arguments into index tuples
     * and to return the number of columns and rows the backend matrix should
     * have.
     */
    template<typename A, typename ... C>
    struct splitter<given, A, C...>
    {
      static constexpr int cols(given&& g, A&& a, C&&... cs)
      {
        return 1;
      }

      static constexpr int rows(given&& g, A&& a, C&&... cs)
      {
        return read_index<A>::read(a) *
        splitter<given,C...>::rows(std::forward<given>(g), std::forward<C>(cs)...);
      }

      static constexpr int posteriors()
      {
      	return 0;
      }

      static constexpr int conditionals()
      {
      	return sizeof...(C)+1;
      }

      typedef vars<> posterior_type;
      typedef vars<A,C...> conditional_type;
      typedef vars<given, A, C...> expanded_type;

      static constexpr bool conditional_distribution = true;

      static typename posterior_type::index_type col_index(const given& g,
						       A&& a, C&&... cs)
      {
      	return std::make_tuple<>();
      }

      static typename conditional_type::index_type row_index(const given& g,
							     A&& a, C&&... cs)
      {
      	return util::tuple::cons(a,splitter<given, C...>::row_index(g, std::forward<C>(cs)...));
      }
    };

		/**
		 * @brief Splits a type list into conditional and posterior variables
		 *
		 * The splitter is the core struct for getting access to the posterior
		 * and conditional variables from a single packed type list (T...). It
		 * also provides methods to convert packed arguments into index tuples
		 * and to return the number of columns and rows the backend matrix should
		 * have.
		 */
    template<typename A>
    struct splitter<given, A>
    {
      static constexpr int cols(given&& g, A&& a)
      {
        return 1;
      }

      static constexpr int rows(given&& g, A&& a)
      {
        return read_index<A>::read(a);
      }

      static constexpr int posteriors()
      {
        return 0;
      }

      static constexpr int conditionals()
      {
        return 1;
      }

      typedef vars<> posterior_type;
      typedef vars<A> conditional_type;
      typedef vars<given, A> expanded_type;

      static constexpr bool conditional_distribution = true;

      static typename posterior_type::index_type col_index(const given& g,
      A&& a)
      {
        return std::make_tuple<>();
      }

      static typename conditional_type::index_type row_index(const given& g,
      A&& a)
      {
        return std::make_tuple(a);
      }
    };

    /**
     * @brief Splits a type list into conditional and posterior variables
     *
     * The splitter is the core struct for getting access to the posterior
     * and conditional variables from a single packed type list (T...). It
     * also provides methods to convert packed arguments into index tuples
     * and to return the number of columns and rows the backend matrix should
     * have.
     */
    template<typename A, typename B, typename D, typename ... C>
    struct splitter<_given<A, B>, D, C...>
    {
      static constexpr int cols(_given<A,B> &&g, D&& d, C&&... cs)
      {
        return read_index<A>::read(g._a);
      }

      static constexpr int rows(_given<A,B> &&g, D&& d, C&&... cs)
      {
        return read_index<B>::read(g._b) *
        splitter<given,D,C...>::rows(prob::given(),
            std::forward<D>(d),
            std::forward<C>(cs)...);
      }

      static constexpr int posteriors()
      {
        return 1;
      }

      static constexpr int conditionals()
      {
        return sizeof...(C)+2;
      }

      typedef vars<A> posterior_type;
      typedef vars<B, D, C...> conditional_type;
      typedef vars<A, given, B, D, C...> expanded_type;

      static constexpr bool conditional_distribution = true;

      static typename posterior_type::index_type col_index(_given<A,B>&& g,
          D&& d, C&&... cs)
      {
        return std::make_tuple(std::forward<A>(g._a));
      }

      static typename conditional_type::index_type row_index(_given<A,B>&& g,
          D&& d, C&&... cs)
      {
        return util::tuple::cons(std::forward<B>(g._b),
            util::tuple::cons(std::forward<D>(d),
                splitter<C...>::row_index(
                    std::forward<C>(cs)...)));
      }

    };

    /**
     * @brief Splits a type list into conditional and posterior variables
     *
     * The splitter is the core struct for getting access to the posterior
     * and conditional variables from a single packed type list (T...). It
     * also provides methods to convert packed arguments into index tuples
     * and to return the number of columns and rows the backend matrix should
     * have.
     */
    template<typename A, typename B>
    struct splitter<_given<A, B>>
    {
      static constexpr int cols(_given<A, B> &&g)
      {
        return read_index<A>::read(g._a);
      }

      static constexpr int rows(_given<A,B> &&g)
      {
        return read_index<B>::read(g._b);
      }

      static constexpr int posteriors()
      {
        return 1;
      }

      static constexpr int conditionals()
      {
        return 1;
      }

      typedef vars<A> posterior_type;
      typedef vars<B> conditional_type;
      typedef vars<A, given, B> expanded_type;

      static constexpr bool conditional_distribution = true;

      static typename posterior_type::index_type col_index(_given<A,B>&& g)
      {
        return std::make_tuple(std::forward<A>(g._a));
      }

      static typename conditional_type::index_type row_index(_given<A,B>&& g)
      {
        return std::make_tuple(std::forward<B>(g._b));
      }

    };

    /**
     * @brief Splits a type list into conditional and posterior variables
     *
     * The splitter is the core struct for getting access to the posterior
     * and conditional variables from a single packed type list (T...). It
     * also provides methods to convert packed arguments into index tuples
     * and to return the number of columns and rows the backend matrix should
     * have.
     */
    template<typename A, typename ... C>
    struct splitter<A, C...>
    {
      static constexpr int cols(A && a, C&&... cs)
      {
        return read_index<A>::read(a) * splitter<C...>::cols(std::forward<C>(cs)...);
      }

      static constexpr int rows(A&& a, C&&... cs)
      {
        return splitter<C...>::rows(std::forward<C>(cs)...);
      }

      static constexpr int posteriors()
      {
        return splitter<C...>::posteriors()+1;
      }

      static constexpr int conditionals()
      {
        return splitter<C...>::conditionals();
      }

      typedef typename util::traits::join<vars, vars<A>,
      typename splitter<C...>::posterior_type>::type posterior_type;
      typedef typename splitter<C...>::conditional_type conditional_type;
      typedef typename util::traits::join<vars, vars<A>,
      typename splitter<C...>::expanded_type>::type expanded_type;

      static constexpr bool conditional_distribution = splitter<C...>::conditional_distribution;

      static typename posterior_type::index_type col_index(A&& a, C&&... cs)
      {
        return util::tuple::cons(std::forward<A>(a),
            splitter<C...>::col_index(std::forward<C>(cs)...));
      }

      static typename conditional_type::index_type row_index(A&& a, C&&... cs)
      {
        return splitter<C...>::row_index(std::forward<C>(cs)...);
      }
    };

    /**
     * @brief Check whether a list of indices is out of bounds
     *
     * Checks whether all indices are within the range of the given
     * number of variables
     *
     * @todo Check for duplicates ...
     */
    template<size_t P, size_t C, size_t Max, size_t HeadI, size_t ...I>
    struct check_indices<P, C, Max, HeadI, I...>
    {
      static constexpr bool valid()
      {
      	// Remaining indices are valid?
        return check_indices<P,C, (HeadI > Max ? HeadI : Max), I...>::valid()
        // This one is within range ?
        		&& HeadI < P+C && HeadI >= 0
        // And if Max > P then check Head > P, no going back from conditionals
        	 	&& (Max >= P ? HeadI >= P : true);
      }
    };

    /**
     * @brief Check whether a list of indices is out of bounds
     *
     * Checks whether all indices are within the range of the given
     * number of variables
     */
    template<size_t P, size_t C, size_t Max, size_t I>
    struct check_indices<P, C, Max, I>
    {
      static constexpr bool valid()
      {
        return I < P + C && I >= 0 && (Max >= P ? I >= P : true);
      }
    };

    /**
     * @brief Adds an offset of one for every index denoting a conditional
     *
     * This struct is needed to take into account that the dummy struct of given
     * takes a place in all index tuples, but is not counting when indexing the
     * variables.
     */
    template<size_t P, int Last, int I>
    struct offset
    {
      static constexpr size_t offset_index()
      {
        //static_assert(I >= Last, "Index selection not strictly monotonic.");
        return I >= P ? I + 1 : I;
      }
    };

    /** @brief Inject a given into a type list on a boolean condition */
    template<typename ...T, template<typename ...> class U, bool Given>
    struct inject_given<U<T...>, Given>
    {
      typedef U<T...> type;
      typedef U<T...> last_type;
    };

    /** @brief Inject a given into a type list on a boolean condition */
    template<typename ...T, template<typename ...> class U>
    struct inject_given<U<T...>, true>
    {
      typedef U<T..., given> type;
      typedef U<given,T...> last_type;
    };

    /** @brief Inject an index into a compile time list on a boolean condition */
    template<size_t ...T, size_t Given, bool GivenCond>
    struct inject_given_index<util::compile_time_list::integer_list<T...>, Given, GivenCond>
    {
      typedef util::compile_time_list::integer_list<T...> type;
    };

    /** @brief Inject an index into a compile time list on a boolean condition */
    template<size_t ...T, size_t Given>
    struct inject_given_index<util::compile_time_list::integer_list<T...>, Given, true>
    {
      typedef util::compile_time_list::integer_list<Given,T...> type;
    };

    /** @brief Workaround to access a probability reference of a distribution via an index tuple
		 *
     */
    template<typename Scalar, typename D, typename T, template<size_t...> class I, size_t... Indices>
    struct tuple_ref_getter<Scalar, D, T, I<Indices...> >
    {
      static Scalar& ref_get(D& d, T&& tuple)
      {
        return d.prob_ref(std::get<Indices>(tuple)...);
      }
    };

    /** @brief Reference access using a subset of the specified arguments
     *
     * The subset is supplied using an index type which is a compile time integer
     * list.
     */
    template<typename Index, typename Scalar, typename D, typename ... Args>
    struct ref
    {
      static Scalar& get(D& d, const Args&... args)
      {
        std::tuple<Args...> arg_tuple = std::make_tuple(args...);
        auto selected_args_tuple = util::tuple::subset(arg_tuple, Index());
        return tuple_ref_getter<
        Scalar, D,
        decltype(selected_args_tuple),
        typename util::compile_time_list::iota_0<Index::size>::type
        >::ref_get(d, std::forward<decltype(selected_args_tuple)>(selected_args_tuple));
      }
    };

    /** @brief Select a subset of types with a supplied index list
     *
     * Select a subset of types with a supplied index list, this directly
     * cares for offseting indices that account for conditional variables
     * and later injecting the \ref given dummy type into the selection.
     *
     * This is not a here be dragons comment, but this is quite heavy type shifting
     * happening here, so if you change stuff, expect it to break quickly and produce
     * dozillions of error messages.
     */
    template<typename ...T, size_t P, size_t C, int Last, int HeadI, template<
        typename ...> class U, int ...I>
    struct indexed_type_selector<U<T...>, P, C, Last, HeadI, I...>
    {
      // Recursively build the type, so  get the tail type
    	typedef typename indexed_type_selector<U<T...>,P,C, HeadI, I...>::result_type tail_type;

    	// Construct the head type by injecting the given if it is needed
      typedef typename inject_given<
      										U<typename util::traits::type_select<
      																	offset<P,Last,HeadI>::offset_index(),
      																	T...>::type>,
      											(Last<(int)P) && (HeadI>=(int)P)>::type head_type;

      // The result type just joins the two type lists
      typedef typename util::traits::join<
          U,
          head_type,
          tail_type>::type result_type;

      // Also create a corresponding index type, that offsets the conditional indices
      typedef typename indexed_type_selector<
      										U<T...>,P,C, HeadI, I...>::index_type::template push_front<
      											offset<P,Last,HeadI>::offset_index()>::type _index_type;

      // And possibly add an index for the given dummy
      typedef typename inject_given_index<_index_type, P,
					  (Last<(int)P) && (HeadI>=(int)P)>::type index_type;

    };

    template<typename ...T, size_t P, size_t C, int Last, int HeadI, template<
        typename ...> class U>
    struct indexed_type_selector<U<T...>, P, C, Last, HeadI>
    {
      typedef typename inject_given<
      		U<typename util::traits::type_select<offset<P,Last,HeadI>::offset_index(), T...>::type>,
				  (Last<(int)P) && (HeadI>=(int)P)>::last_type result_type;

      typedef typename inject_given_index<
          util::compile_time_list::integer_list<offset<P,Last,HeadI>::offset_index()>,
          P,
          (Last<(int)P) && (HeadI>=(int)P)>::type index_type;
    };


    template<size_t ...A, typename B, size_t P, size_t Head>
    struct cased_append<util::compile_time_list::integer_list<A...>, B, P, true, Head>
    {
      typedef util::compile_time_list::integer_list<Head, A...> col_index_type;
      typedef B row_index_type;
    };

    template<typename A, size_t ...B, size_t P, size_t Head>
    struct cased_append<A, util::compile_time_list::integer_list<B...>, P, false, Head>
    {
      typedef A col_index_type;
      typedef util::compile_time_list::integer_list<Head-P, B...> row_index_type;
    };

    template<typename T, size_t P, size_t C, int HeadI, int ...I>
    struct index_splitter<T, P, C, HeadI, I...>
    {
      typedef index_splitter<T,P,C, I...> tail_struct;
      typedef cased_append<typename tail_struct::col_index_type,
      					typename tail_struct::row_index_type,
      					P,
      					(HeadI<P), HeadI> append_type;
      typedef typename append_type::col_index_type col_index_type;
      typedef typename append_type::row_index_type row_index_type;
    };

    template<typename T, size_t P, size_t C, int HeadI>
    struct index_splitter<T, P, C, HeadI>
    {
      typedef cased_append<util::compile_time_list::integer_list<>,
          util::compile_time_list::integer_list<>, P, (HeadI < P), HeadI> append_type;
      typedef typename append_type::col_index_type col_index_type;
      typedef typename append_type::row_index_type row_index_type;

    };

  }
}

#endif /* _SPLITTER_H_ */
