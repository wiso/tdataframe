#ifndef TTMPDATAFRAME
#define TTMPDATAFRAME

#include "TDataFrameTypes.hpp"
#include "helpers.hpp"
#include "metautils.hpp"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include <string>
#include <tuple>
#include <list>
#include <memory> //std::static_pointer_cast


template<typename Filter, typename PrevData>
class TTmpDataFrame {
   template<typename A, typename B> friend class TTmpDataFrame;
   using filter_types = typename f_traits<Filter>::arg_types_tuple;
   using filter_ind = typename gens<std::tuple_size<filter_types>::value>::type;

   public:
   TTmpDataFrame(const BranchList& _bl, Filter _f, PrevData& _pd)
      : bl(_bl), f(_f), pd(_pd), def_bl(_pd.def_bl), dir(_pd.dir),
        tree_name(_pd.tree_name) {}

   template<typename NewFilter>
   auto filter(NewFilter f, const BranchList& bl = {})
      -> TTmpDataFrame<NewFilter, decltype(*this)>;

   template<typename F>
   void foreach(const BranchList& branches, F f);

   template<typename T>
   TH1F fillhist(std::string branch, unsigned nbins = 100, std::string suffix = "");

   std::list<unsigned> collect_entries();

   template<typename T>
   std::list<T> get(std::string branch);

   private:
   bool apply_filters() {
      return apply_filters(filter_types(), filter_ind());
   }

   template<int... S, typename... types>
   bool apply_filters(std::tuple<types...>, seq<S...>);

   template<typename F, int... S, typename... types>
   void loop_and_apply(TTreeReader& r, F f, const BranchList& branches,
                       std::tuple<types...> types_tuple, seq<S...> intseq);

   template<typename F, int... S, typename... types>
   void foreach_helper(const BranchList& branches, F f,
                       std::tuple<types...> types_tuple, seq<S...> intseq);

   void build_filter_tvb(TTreeReader& r);

   const BranchList& bl;
   Filter f;
   PrevData pd;
   const BranchList& def_bl;
   TDirectory* dir;
   const std::string& tree_name;
   TVBVec filter_tvb;
};


template<typename Filter, typename PrevData>
template<typename NewFilter>
auto TTmpDataFrame<Filter,PrevData>::filter(NewFilter f, const BranchList& bl)
-> TTmpDataFrame<NewFilter, decltype(*this)>
{
   bool use_def_bl = check_filter(f, bl, def_bl);
   const BranchList& actual_bl = use_def_bl ? def_bl : bl;
   // Create a TTmpDataFrame that contains *this (and the new filter)
   return TTmpDataFrame<NewFilter, decltype(*this)>(actual_bl, f, *this);
}


template<typename Filter, typename PrevData>
template<typename F>
void TTmpDataFrame<Filter,PrevData>::foreach(const BranchList& branches, F f)
{
   using arg_t_tuple = typename f_traits<F>::arg_types_tuple;
   using arg_indexes = typename gens<std::tuple_size<arg_t_tuple>::value>::type;
   #ifdef R__USE_IMT
      if(ROOT::IsImplicitMTEnabled()) {
         auto file_name = dir->GetName();
         ROOT::TTreeProcessor tp(file_name, tree_name);
         tp.Process(
            [&f, &branches, this] (TTreeReader& my_r) {
               this->loop_and_apply(my_r, f, branches, arg_t_tuple(), arg_indexes());
            }
         );
      } else {
         TTreeReader r(tree_name.c_str(), dir);
         if(!check_reader(r))
            return;
         loop_and_apply(r, f, branches, f_arg_types(), f_arg_indexes());
      }
   #else
      TTreeReader r(tree_name.c_str(), dir);
      if(!check_reader(r))
         return;
      loop_and_apply(r, f, branches, arg_t_tuple(), arg_indexes());
   #endif
}


template<typename Filter, typename PrevData>
std::list<unsigned> TTmpDataFrame<Filter,PrevData>::collect_entries()
{
   std::list<unsigned> l;
   TTreeReader r(tree_name.c_str(), dir);
   if(check_reader(r)) {
      build_filter_tvb(r);
      while(r.Next())
         if(apply_filters())
            l.push_back(r.GetCurrentEntry());
   }

   return l;
}

template<typename Filter, typename PrevData>
template<typename T>
std::list<T> TTmpDataFrame<Filter,PrevData>::get(std::string branch)
{
   TTreeReader r(tree_name.c_str(), dir);
   std::list<T> res;
   if(check_reader(r)) {
      build_filter_tvb(r);
      TTreeReaderValue<T> v(r, branch.c_str());
      while(r.Next())
         if(apply_filters())
            res.push_back(*v);
   }

   return res;
}


template<typename Filter, typename PrevData>
template<typename T>
TH1F TTmpDataFrame<Filter,PrevData>::fillhist(std::string branch, unsigned nbins, std::string suffix)
{
   // histogram with automatic binning
   TH1F h(("h_" + branch + suffix).c_str(), branch.c_str(), nbins, 0., 0.);
   TTreeReader r(tree_name.c_str(), dir);
   if(check_reader(r)) {
      build_filter_tvb(r);
      TTreeReaderValue<T> v(r, branch.c_str());
      while(r.Next())
         if(apply_filters())
            h.Fill(*v);
   }

   return h;
}


template<typename Filter, typename PrevData>
template<int... S, typename... types>
bool TTmpDataFrame<Filter,PrevData>::apply_filters(std::tuple<types...>, seq<S...>)
{
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


template<typename Filter, typename PrevData>
template<typename F, int... S, typename... types>
void TTmpDataFrame<Filter,PrevData>::loop_and_apply(
   TTreeReader& r, F f, const BranchList& branches,
   std::tuple<types...> types_tuple, seq<S...> intseq)
{
   auto f_tvb = build_tvb(r, branches, types_tuple, intseq);
   build_filter_tvb(r);
   while(r.Next()) {
      if(apply_filters())
         f(*(std::static_pointer_cast<TTreeReaderValue<types>>(f_tvb[S]))->Get()...);
   }
}


template<typename Filter, typename PrevData>
void TTmpDataFrame<Filter,PrevData>::build_filter_tvb(TTreeReader& r) {
   filter_tvb = build_tvb(r, bl, filter_types(), filter_ind());
   pd.build_filter_tvb(r);
   return;
}

#endif // TTMPDATAFRAME
