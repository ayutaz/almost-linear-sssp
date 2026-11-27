#include "paper/model/random.hpp"

namespace sssp::paper::model {

// Instantiations for the library generator (the Sec 6 sampling code uses these).
template std::uint64_t uniform_below<WordRng>(WordRng&, std::uint64_t);
template std::vector<std::uint32_t> random_permutation<WordRng>(std::size_t, WordRng&);
template void random_shuffle<std::uint32_t, WordRng>(std::span<std::uint32_t>, WordRng&);

}  // namespace sssp::paper::model
