#pragma once

#include <string_view>

namespace sssp {

/// Library semantic version.
std::string_view version_string() noexcept;

/// arXiv identifier (with version) of the paper this project reproduces.
std::string_view paper_arxiv_id() noexcept;

}  // namespace sssp
