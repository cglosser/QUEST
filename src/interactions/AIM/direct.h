#pragma once

#include "interactions/AIM/grid.h"
#include "interactions/history_interaction.h"

namespace AIM {
  class DirectInteraction;
}

class AIM::DirectInteraction final : public HistoryInteraction {
 public:
  DirectInteraction(
      std::shared_ptr<const DotVector>,
      std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
      Propagation::Kernel<cmplx> &,
      const int,
      const double,
      const double,
      std::shared_ptr<const std::vector<Grid::ipair_t>>);

  const ResultArray &evaluate(const int) final;
  const ResultArray &evaluate_present_field(const int) final;
  boost::multi_array<cmplx, 2> coefficient_table(Propagation::Kernel<cmplx> &,
                                                 std::vector<int> &) const;

 private:
  std::shared_ptr<const std::vector<Grid::ipair_t>> interaction_pairs_;
  std::array<int, 2> shape_;
  std::vector<int> floor_delays_;
  boost::multi_array<cmplx, 2> coefficients_;
};
