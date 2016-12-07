#ifndef TDATAFRAME
#define TDATAFRAME

#include <vector>
#include <list>
#include <tuple>
#include <functional>
#include <string>
#include <type_traits> //std::is_same
#include <stdexcept> //std::runtime_error
#include <memory> //std::shared_ptr
#include "TTree.h"
#include "TH1F.h"
#include "TTreeReaderValue.h"
#include "TTreeReader.h"
#include "ROOT/TTreeProcessor.hxx"
#include "TFile.h"

/******* meta-utils **********/
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
// credits to Johannes Schaub (see http://stackoverflow.com/q/7858817/2195972)
template<int ...>
struct seq {};

template<int N, int ...S>
struct gens : gens<N-1, N-1, S...> {};

template<int ...S>
struct gens<0, S...>{
   typedef seq<S...> type;
};
/****************************/


using BranchList = std::vector<std::string>;
using EntryList = std::list<unsigned>;
using TVBVec = std::vector<std::shared_ptr<ROOT::Internal::TTreeReaderValueBase>>;


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


// forward declaration (needed by TDataFrame)
template<class Filter, class PrevData>
class TTmpDataFrame;


class TDataFrame {
   template<class A, class B> friend class TTmpDataFrame;

   public:
   TDataFrame(TTreeReader& _t, const BranchList& _bl = {}) : t(_t), def_bl(_bl) {}

   template<class Filter>
   auto filter(Filter f, const BranchList& bl = {}) -> TTmpDataFrame<Filter, decltype(*this)> {
      bool use_def_bl = check_filter(f, bl, def_bl);
      // Every time this TDataFrame is (re)used we want a fresh TTreeReader
      t.Restart();
      // Create a TTmpDataFrame that contains *this (and the new filter)
      const BranchList& actual_bl = use_def_bl ? def_bl : bl;
      return TTmpDataFrame<Filter, decltype(*this)>(t, actual_bl, f, *this);
   }

   private:
   bool apply_filters() {
      // Dummy call (end of recursive chain of calls)
      return true;
   }

   void build_filter_tvb(TTreeReader&) {
      // Dummy call (end of recursive chain of calls)
      return;
   }

   TTreeReader& t;
   //! the branchlist to fall back to if none is specified in filters and actions
   const BranchList def_bl;
};


template<int... S, typename... arg_types>
TVBVec build_tvb(TTreeReader& t, const BranchList& bl, std::tuple<arg_types...>, seq<S...>) {
   // Build vector of pointers to TTreeReaderValueBase.
   // tvb[i] is a TTreeReaderValue specialized for the i-th arg_type
   // S must be a sequence of sizeof...(arg_types) integers
   // arg_types and S are expanded simultaneously by "..."
   TVBVec tvb{ std::make_shared<TTreeReaderValue<arg_types>>(t, bl.at(S).c_str())... };
   return tvb;
}


template<class Filter, class PrevData>
class TTmpDataFrame {
   template<class A, class B> friend class TTmpDataFrame;
   template<class A, class B> friend class TTmpDataProcessor;
   using filter_types = typename f_traits<Filter>::arg_types_tuple;
   using filter_ind = typename gens<std::tuple_size<filter_types>::value>::type;

   public:
   TTmpDataFrame(TTreeReader& _t, const BranchList& _bl, Filter _f, PrevData& _pd)
      : t(_t), bl(_bl), f(_f), pd(_pd), def_bl(_pd.def_bl) {}

   template<class NewFilter>
   auto filter(NewFilter f, const BranchList& bl = {}) -> TTmpDataFrame<NewFilter, decltype(*this)> {
      bool use_def_bl = check_filter(f, bl, def_bl);
      const BranchList& actual_bl = use_def_bl ? def_bl : bl;

      // Create a TTmpDataFrame that contains *this (and the new filter)
      return TTmpDataFrame<NewFilter, decltype(*this)>(t, actual_bl, f, *this);
   }

   EntryList collect_entries() {
      EntryList l;
      build_filter_tvb(t);
      while(t.Next())
         if(apply_filters())
            l.push_back(t.GetCurrentEntry());
      return l;
   }

   template<class T>
   std::list<T> get(std::string branch) {
      std::list<T> res;
      build_filter_tvb(t);
      TTreeReaderValue<T> v(t, branch.c_str());
      while(t.Next())
         if(apply_filters())
            res.push_back(*v);

      return res;
   }

   template<class F>
   void foreach(const BranchList& branches, F f) {
      using f_arg_types = typename f_traits<F>::arg_types_tuple;
      using f_arg_indexes = typename gens<std::tuple_size<f_arg_types>::value>::type;
      apply_function(branches, f, f_arg_types(), f_arg_indexes());
   }

   template<class T>
   TH1F fillhist(std::string branch, unsigned nbins = 100,
                 std::string name_suffix = "") {
      // histogram with automatic binning
      TH1F h(("h_" + branch).c_str(), branch.c_str(), nbins, 0., 0.);
      build_filter_tvb(t);
      TTreeReaderValue<T> v(t, branch.c_str());
      while(t.Next())
         if(apply_filters())
            h.Fill(*v);

      return h;
   }

   private:
   bool apply_filters() {
      return apply_filters(filter_types(), filter_ind());
   }

   template<int... S, typename... types>
   bool apply_filters(std::tuple<types...>, seq<S...>) {
      // Recursive call to all previous filters
      if(!pd.apply_filters())
         return false;

      // Apply our filter
      // Take each pointer in tvb, cast it to a pointer to the
      // correct specialization of TTreeReaderValue, and get its content.
      // S expands to a sequence of integers 0 to sizeof...(types)-1
      // S and types are expanded simultaneously by "..."
      return f(*(std::static_pointer_cast<TTreeReaderValue<types>>(filter_tvb[S]))->Get() ...);
   }

   template<typename F, int... S, typename... types>
   void loop_and_apply(TTreeReader& my_t, F f, const BranchList& branches,
                       std::tuple<types...> types_tuple, seq<S...> intseq) {
      auto f_tvb = build_tvb(my_t, branches, types_tuple, intseq);
      build_filter_tvb(my_t); // recursive call to build TTree
      while(my_t.Next()) {
         if(apply_filters())
            f(*(std::static_pointer_cast<TTreeReaderValue<types>>(f_tvb[S]))->Get()...);
      }
   }

   template<typename F, int... S, typename... types>
   // TODO this should probably be called foreach_helper
   void apply_function(const BranchList& branches, F f, std::tuple<types...> types_tuple, seq<S...> intseq) {
#ifdef R__USE_IMT
      if(ROOT::IsImplicitMTEnabled()) {
         auto tree = t.GetTree();
         auto tree_name = tree->GetName();
         auto file = tree->GetCurrentFile();
         auto file_name = file->GetName();
         ROOT::TTreeProcessor tp(file_name, tree_name);
         tp.Process([&f, &branches, &types_tuple, &intseq, this] (TTreeReader& my_t) -> void {
            build_filter_tvb(my_t);
            this->loop_and_apply(my_t, f, branches, types_tuple, intseq);
         });
      } else {
         build_filter_tvb(t);
         loop_and_apply(t, f, branches, types_tuple, intseq);
      }
#else
      build_filter_tvb(t);
      loop_and_apply(t, f, branches, types_tuple, intseq);
#endif // R__USE_IMT
   }

   void build_filter_tvb(TTreeReader& my_t) {
      filter_tvb = build_tvb(my_t, bl, filter_types(), filter_ind());
      pd.build_filter_tvb(my_t);
      return;
   }

   TTreeReader& t;
   const BranchList& bl;
   Filter f;
   PrevData pd;
   const BranchList& def_bl;
   TVBVec filter_tvb;
};

#endif // TDATAFRAME
