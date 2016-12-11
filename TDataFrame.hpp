#ifndef TDATAFRAME
#define TDATAFRAME

#include "helpers.hpp" // check_filter
#include "TDataFrameTypes.hpp" // BranchList, EntryList, TVBVec
#include "TTmpDataFrame.hpp"
#include "TDirectory.h"
#include <string>


class TDataFrame {
   template<typename A, typename B> friend class TTmpDataFrame;

   public:
   TDataFrame(const std::string _tree_name, TDirectory* _dir = nullptr, const BranchList& _bl = {})
      : tree_name(_tree_name), dir(_dir), def_bl(_bl)
   {
      if (dir == nullptr)
         dir = gDirectory;
   }

   template<typename Filter>
   auto filter(Filter f, const BranchList& bl = {})
      -> TTmpDataFrame<Filter, decltype(*this)>;

   private:
   //! Dummy call (end of recursive chain of calls)
   bool apply_filters() { return true; }

   //! Dummy call (end of recursive chain of calls)
   void build_filter_tvb(TTreeReader&) { return; }

   //! Name of the TTree to process
   const std::string tree_name;
   //! Directory (e.g. TFile) where the TTree is located. Defaults to current directory if null.
   TDirectory* dir;
   //! BranchList to fall back to if none is specified in filters/actions
   const BranchList def_bl;
};


template<typename Filter>
auto TDataFrame::filter(Filter f, const BranchList& bl)
-> TTmpDataFrame<Filter, decltype(*this)>
{
   bool use_def_bl = check_filter(f, bl, def_bl);
   const BranchList& actual_bl = use_def_bl ? def_bl : bl;
   // Create a TTmpDataFrame that contains *this (and the new filter)
   return TTmpDataFrame<Filter, decltype(*this)>(actual_bl, f, *this);
}

#endif // TDATAFRAME
