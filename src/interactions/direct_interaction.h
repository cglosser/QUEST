#ifndef DIRECT_INTERACTION_H
#define DIRECT_INTERACTION_H

#include "history_interaction.h"

class DirectInteraction final : public HistoryInteraction {
 public:
  DirectInteraction(
      std::shared_ptr<const DotVector>,
      std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
      Propagation::Kernel<cmplx> &,
      const int,
      const double,
      const double);

  const ResultArray &evaluate(const int) final;
  const ResultArray &evaluate_present_field(const int) final;

 private:
  int num_interactions;
  std::vector<int> floor_delays;
  boost::multi_array<cmplx, 2> coefficients;

  void build_coefficient_table(Propagation::Kernel<cmplx> &);

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
};

#endif
