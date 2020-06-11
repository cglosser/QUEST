#ifndef FARFIELD_H
#define FARFIELD_H

#include "aim_base.h"
#include "fourier.h"

namespace AIM {
  class Farfield;
}

class AIM::Farfield final : public AimBase {
 public:
  Farfield(const std::shared_ptr<const DotVector>,
           const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
           const int,
           const double,
           const double,
           std::shared_ptr<const Grid>,
           std::shared_ptr<const Expansions::ExpansionTable>,
           Expansions::ExpansionFunction,
           Normalization::SpatialNorm);
  ~Farfield() = default;

  const ResultArray &evaluate_present_field(const int step) final
  {
    fill_source_table(step);
    propagate_present_field(step);
    fill_results_table(step);

    return results;
  }
  const ResultArray &evaluate(const int step) final
  {
    fill_source_table(step);
    propagate(step);
    fill_results_table(step);
    return results;
  }

 private:
  std::array<int, 4> table_dimensions_;

  // The propagation table corresponds to delta(t - R/c)/R and thus holds
  // *scalar* quantities; the source and obs tables correspond to currents and
  // fields and thus hold *vector* quantities
  spacetime::vector<cmplx> propagation_table_;
  spacetime::vector3d<cmplx> source_table_, obs_table_;

  Eigen::Array3Xcd temp_observers;

  TransformPair spatial_vector_transforms_;

  spacetime::vector<cmplx> make_propagation_table() const;
  void fill_source_table(const int);
  void propagate(const int);
  void propagate_present_field(const int);
  void fill_results_table(const int);

  void fill_gmatrix_table(spacetime::vector<cmplx> &) const;
  TransformPair spatial_fft_plans();
};

#endif
