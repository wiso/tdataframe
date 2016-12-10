#ifndef HELPERS
#define HELPERS

#include <tuple>
#include <string> // std::to_string
#include <type_traits> //std::is_same
#include <memory> //std::make_shared
#include "metautils.hpp" // f_traits, seq
#include "TDataFrameTypes.hpp" //BranchList
#include "TTreeReader.h"

// Build vector of pointers to TTreeReaderValueBase.
// tvb[i] is a TTreeReaderValue specialized for the i-th arg_type
// S must be a sequence of sizeof...(arg_types) integers
// arg_types and S are expanded simultaneously by "..."
template<int... S, typename... arg_types>
TVBVec build_tvb(TTreeReader& t, const BranchList& bl, std::tuple<arg_types...>, seq<S...>) {
   TVBVec tvb{ std::make_shared<TTreeReaderValue<arg_types>>(t, bl.at(S).c_str())... };
   return tvb;
}


// Check that filter is valid and is a good match for the branch list.
// return true if the default branch list should be used for this filter,
// false otherwise.
template<typename Filter>
bool check_filter(Filter f, const BranchList& bl, const BranchList& def_bl) {
   using filter_ret_t = typename f_traits<Filter>::ret_t;
   static_assert(std::is_same<filter_ret_t, bool>::value,
                 "filter functions must return a bool");
   using filter_args_tuple = typename f_traits<Filter>::arg_types_tuple;
   auto n_args = std::tuple_size<filter_args_tuple>::value;

   bool use_def_bl = false;
   if(n_args != bl.size()) {
      if(bl.size() == 0 && n_args == def_bl.size()) {
         use_def_bl = true;
      } else {
         auto msg = "mismatch between number of filter arguments (" \
                     + std::to_string(n_args) + ") and number of branches (" \
                     + std::to_string(bl.size()?:def_bl.size()) + ")";
         throw std::runtime_error(msg);
      }
   }

   return use_def_bl;
}

#endif // HELPERS
