#ifndef METAUTILS
#define METAUTILS

#include <tuple>

// extract parameter types from a callable object
template<typename T>
struct f_traits : public f_traits<decltype(&T::operator())> {};

// lambdas and std::function
template<typename R, typename T, typename... Args>
struct f_traits<R(T::*)(Args...) const> {
   using arg_types_tuple = typename std::tuple<Args...>;
   using ret_t = R;
};

// mutable lambdas and functor classes
template<typename R, typename T, typename... Args>
struct f_traits<R(T::*)(Args...)> {
   using arg_types_tuple = typename std::tuple<Args...>;
   using ret_t = R;
};

// free functions
template<typename R, typename... Args>
struct f_traits<R(*)(Args...)> {
   using arg_types_tuple = typename std::tuple<Args...>;
   using ret_t = R;
};


// compile-time integer sequence generator
// e.g. calling gens<3>::type() instantiates a seq<0,1,2>
template<int ...>
struct seq {};

template<int N, int ...S>
struct gens : gens<N-1, N-1, S...> {};

template<int ...S>
struct gens<0, S...>{
   typedef seq<S...> type;
};

#endif // METAUTILS
