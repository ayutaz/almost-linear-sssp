// Explicit instantiations of the split-findmin structures (Sec 3.5) for the two key types used
// by the evaluator: exact integer lengths and double lengths.

#include "paper/primitives/split_findmin.hpp"

#include <cstdint>

namespace sssp::paper {

template class SplitFindminRef<std::uint64_t>;
template class SplitFindminRef<double>;
template class DyadicSplitFindmin<std::uint64_t>;
template class DyadicSplitFindmin<double>;
template class BlockedSplitFindmin<DyadicSplitFindmin<std::uint64_t>,
                                   DyadicSplitFindmin<std::uint64_t>>;
template class BlockedSplitFindmin<DyadicSplitFindmin<double>, DyadicSplitFindmin<double>>;
template class BlockedSplitFindmin<DyadicSplitFindmin<std::uint64_t>,
                                   OneLevelBlockedSplitFindmin<std::uint64_t>>;
template class BlockedSplitFindmin<DyadicSplitFindmin<double>, OneLevelBlockedSplitFindmin<double>>;

}  // namespace sssp::paper
