#ifndef NEARFIELD_H
#define NEARFIELD_H

#include <limits>
#include "aim_base.h"

namespace AIM {
  class Nearfield;
}

class AIM::Nearfield final : public AimBase {
 public:
  Nearfield(const std::shared_ptr<const DotVector>,
            const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
            const int,
            const double,
            const double,
            std::shared_ptr<const Grid>,
            std::shared_ptr<const Expansions::ExpansionTable>,
            Expansions::ExpansionFunction,
            Normalization::SpatialNorm,
            std::shared_ptr<const std::vector<Grid::ipair_t>>,
            const double = 0);
  ~Nearfield() = default;

  const ResultArray &evaluate(const int) final;
  const ResultArray &evaluate_present_field(const int) final;

 private:
  struct support_range_t {
    int begin, end;
  };

  double omega_;
  std::shared_ptr<const std::vector<Grid::ipair_t>> interaction_pairs_;
  std::array<int, 3> shape_;
  std::vector<support_range_t> support_;
  boost::multi_array<cmplx, 3> coefficients_;

  boost::multi_array<cmplx, 3> coefficient_table();
};

#endif
