#include "paper/model/length.hpp"

namespace sssp::paper::model {

// The two weight types used throughout the project (uint64_t exact, double IEEE-754).
template class Len<std::uint64_t>;
template class Len<double>;
template class SignedPair<std::uint64_t>;
template class SignedPair<double>;
template class LenOps<std::uint64_t>;
template class LenOps<double>;

}  // namespace sssp::paper::model
