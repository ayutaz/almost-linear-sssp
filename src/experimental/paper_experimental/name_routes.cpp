// Mode B entry to NameRoutes (arXiv 2611.34911v1, Sec 6.1 Alg 5, Sec 6.2 l.1491), task
// P5-EXP-NAMEROUTES. See name_routes.hpp for the contract.
#include "experimental/paper_experimental/name_routes.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>

#include "paper/model/params.hpp"
#include "paper/pivots/estimate.hpp"

namespace sssp::paper_experimental {
namespace {

constexpr std::uint64_t kWordMax = ~std::uint64_t{0};
constexpr std::string_view kCapPrefix = "name_routes: retry cap of ";

}  // namespace

bool NameRoutesSetting::is_paper() const {
  return p_log2 == paper_p_log2 && h_log2 == paper_h_log2 && sample_size == paper_sample_size &&
         residual_limit == paper_residual_limit && !retry_cap.has_value();
}

NameRoutesSetting name_routes_setting(const PaperParams& params, std::uint64_t t) {
  if (t >= (std::uint64_t{1} << 32)) throw std::length_error("name_routes_setting: t >= 2^32");
  const OuterParams paper_outer;
  // For t < 2 there are no edges and Estimate is never called (l.1422); H is reported at t = 1.
  const std::uint64_t th = std::max<std::uint64_t>(t, 1);
  const paper::pivots::NameRoutesParams used = params.name_routes_params(th);
  const auto t_size = static_cast<std::size_t>(t);
  NameRoutesSetting s;
  s.t = t;
  s.paper_p_log2 = paper::pivots::name_routes_constants(t_size).p_log2;  // Alg 5 l.2-3
  s.p_log2 = used.p_log2.value_or(s.paper_p_log2);                       // or D-17
  const auto ec = paper::pivots::estimate_constants(static_cast<std::size_t>(th));
  s.paper_h_log2 = ec.eta;                           // Sec 6.2 l.1491
  s.h_log2 = used.estimate.h_log2.value_or(ec.eta);  // or D-9
  // Alg 5 l.5-6 with the products capped as in the shared component.
  s.sample_size = paper::model::capped_mul(used.sample_factor, t, kWordMax);
  s.paper_sample_size = paper::model::capped_mul(paper_outer.sample_factor, t, kWordMax);
  s.residual_limit = paper::model::capped_mul(used.residual_factor, t, kWordMax);
  s.paper_residual_limit = paper::model::capped_mul(paper_outer.residual_factor, t, kWordMax);
  s.retry_cap = used.retry_cap;  // D-10
  return s;
}

template <class W>
ExperimentalRouteNames<W> name_routes(const paper::pivots::PivotRoutingInstance<W>& instance,
                                      paper::model::WordRng& rng, const PaperParams& params) {
  const std::uint64_t t = instance.t;
  ExperimentalRouteNames<W> out;
  out.setting = name_routes_setting(params, t);
  // Alg 5 itself: the shared Mode A component with the derived parameters (D-9 needs t >= 1;
  // see name_routes_setting for t = 0).
  out.names = paper::pivots::name_routes(instance, rng,
                                         params.name_routes_params(std::max<std::uint64_t>(t, 1)));
  if (out.names.p_log2 != out.setting.p_log2) {
    throw std::logic_error("paper_experimental::name_routes: P differs from the setting");
  }
  // D-10: when the cap stops the loop the shared component counts "name_routes.retry_cap_hits"
  // (in every build) and pushes one FallbackRecord whose edge count is the residual that
  // certified its own weight. A hit without its record would be a silent fallback.
  const auto hits = out.names.diag.counters.find("name_routes.retry_cap_hits");
  out.retry_cap_reached = hits != out.names.diag.counters.end() && hits->second > 0;
  std::uint64_t records = 0;
  for (const FallbackRecord& f : out.names.diag.fallbacks) {
    if (f.reason.starts_with(kCapPrefix)) {
      ++records;
      out.fallback_edges += f.edges;
    }
  }
  if (records != (out.retry_cap_reached ? 1U : 0U)) {
    throw std::logic_error("paper_experimental::name_routes: retry cap hit and record disagree");
  }
  return out;
}

template ExperimentalRouteNames<std::uint64_t> name_routes<std::uint64_t>(
    const paper::pivots::PivotRoutingInstance<std::uint64_t>&, paper::model::WordRng&,
    const PaperParams&);
template ExperimentalRouteNames<double> name_routes<double>(
    const paper::pivots::PivotRoutingInstance<double>&, paper::model::WordRng&, const PaperParams&);

}  // namespace sssp::paper_experimental
