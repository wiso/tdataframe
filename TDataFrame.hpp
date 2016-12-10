#ifndef TDATAFRAME
#define TDATAFRAME

#include "helpers.hpp" // check_filter
#include "TDataFrameTypes.hpp" // BranchList, EntryList, TVBVec
#include "TTmpDataFrame.hpp"
#include "TTreeReader.h"


class TDataFrame {
   template<typename A, typename B> friend class TTmpDataFrame;

   public:
   TDataFrame(TTreeReader& _t, const BranchList& _bl = {})
      : t(_t), def_bl(_bl) {}

   template<typename Filter>
   auto filter(Filter f, const BranchList& bl = {})
      -> TTmpDataFrame<Filter, decltype(*this)>;

   private:
   //! Dummy call (end of recursive chain of calls)
   bool apply_filters() { return true; }

   //! Dummy call (end of recursive chain of calls)
   void build_filter_tvb(TTreeReader&) { return; }

   TTreeReader& t;
   //! the branchlist to fall back to if none is specified
   const BranchList def_bl;
};


template<typename Filter>
auto TDataFrame::filter(Filter f, const BranchList& bl)
-> TTmpDataFrame<Filter, decltype(*this)>
{
   // Every time this TDataFrame is (re)used we want a fresh TTreeReader
   t.Restart();
   bool use_def_bl = check_filter(f, bl, def_bl);
   const BranchList& actual_bl = use_def_bl ? def_bl : bl;
   // Create a TTmpDataFrame that contains *this (and the new filter)
   return TTmpDataFrame<Filter, decltype(*this)>(t, actual_bl, f, *this);
}

#endif // TDATAFRAME
