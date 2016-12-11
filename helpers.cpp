#include "helpers.hpp"
#include "TError.h"

bool check_reader(const TTreeReader& r) {
   if(r.IsZombie()) {
      Error("check_reader", "could not build TTreeReader from desired tree, aborting action");
      return false;
   }
   
   return true;
}
