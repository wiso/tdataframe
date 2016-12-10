#ifndef TDATAFRAMETYPES
#define TDATAFRAMETYPES

#include <vector>
#include <string>
#include <memory>
#include "TTreeReaderValue.h"

using BranchList = std::vector<std::string>;
using TVBVec = std::vector<std::shared_ptr<ROOT::Internal::TTreeReaderValueBase>>;

#endif // TDATAFRAMETYPES
