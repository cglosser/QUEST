#include "nearfield.h"

AIM::Nearfield::Nearfield(
    const std::shared_ptr<const DotVector> dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    const int interp_order,
    const double c0,
    const double dt,
    std::shared_ptr<const Grid> grid,
    std::shared_ptr<const Expansions::ExpansionTable> expansion_table,
    Expansions::ExpansionFunction expansion_function,
    Normalization::SpatialNorm normalization,
    std::shared_ptr<const std::vector<Grid::ipair_t>> interaction_pairs,
    const double omega)
    : AimBase(dots,
              history,
              interp_order,
              c0,
              dt,
              grid,
              expansion_table,
              expansion_function,
              normalization),
      omega_{omega},
      interaction_pairs_{std::move(interaction_pairs)},
      shape_({{static_cast<int>(interaction_pairs_->size()),
               grid->max_transit_steps(c0, dt) + 12, 2}}),
      support_(shape_[0], {std::numeric_limits<int>::max(), 0}),
      coefficients_{coefficient_table()}
{
}

const InteractionBase::ResultArray &AIM::Nearfield::evaluate(const int time_idx)
{
  constexpr int RH0_01 = 1;

  results.setZero();
  past_terms_of_convolution.setZero();

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];

    for(int t = support_[pair_idx].begin + 1; t < support_[pair_idx].end; ++t) {
      const int s = std::max(
          time_idx - t, static_cast<int>(history->array_.index_bases()[1]));

      past_terms_of_convolution[pair.first] +=
          (history->array_[pair.second][s][0])[RHO_01] *
          coefficients_[pair_idx][t][0];
      past_terms_of_convolution[pair.second] +=
          (history->array_[pair.first][s][0])[RHO_01] *
          coefficients_[pair_idx][t][1];
    }
    int t = support_[pair_idx].begin;
    const int s = std::max(time_idx - t,
                           static_cast<int>(history->array_.index_bases()[1]));
    assert(support_[pair_idx].begin == 0);

    results[pair.first] += (history->array_[pair.second][s][0])[RHO_01] *
                           coefficients_[pair_idx][t][0];
    results[pair.second] += (history->array_[pair.first][s][0])[RHO_01] *
                            coefficients_[pair_idx][t][1];
  }
  results += past_terms_of_convolution;

  return results;
}

const InteractionBase::ResultArray &AIM::Nearfield::evaluate_present_field(
    const int time_idx)
{
  constexpr int RH0_01 = 1;

  results.setZero();

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];

    int t = support_[pair_idx].begin;
    const int s = std::max(time_idx - t,
                           static_cast<int>(history->array_.index_bases()[1]));
    assert(support_[pair_idx].begin == 0);

    results[pair.first] += (history->array_[pair.second][s][0])[RHO_01] *
                           coefficients_[pair_idx][t][0];
    results[pair.second] += (history->array_[pair.first][s][0])[RHO_01] *
                            coefficients_[pair_idx][t][1];
  }
  results += past_terms_of_convolution;

  return results;
}

boost::multi_array<cmplx, 3> AIM::Nearfield::coefficient_table()
{
  boost::multi_array<cmplx, 3> coefficients(shape_);
  std::fill(coefficients.data(),
            coefficients.data() + coefficients.num_elements(), cmplx(0, 0));

  Interpolation::DerivFive lagrange(dt);

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];
    const auto &dot0 = (*dots)[pair.first], &dot1 = (*dots)[pair.second];

    for(size_t i = 0; i < expansion_table->shape()[1]; ++i) {
      const auto &e0 = (*expansion_table)[pair.first][i];

      for(size_t j = 0; j < expansion_table->shape()[1]; ++j) {
        const auto &e1 = (*expansion_table)[pair.second][j];

        if(e0.index == e1.index) continue;

        Eigen::Vector3d dr(grid->spatial_coord_of_box(e1.index) -
                           grid->spatial_coord_of_box(e0.index));

        // == Retardation quantities ==

        const double arg = dr.norm() / (c0 * dt);
        const auto split_arg = split_double(arg);

        support_[pair_idx].begin =
            std::min(support_[pair_idx].begin, split_arg.first);
        support_[pair_idx].end = std::max(
            support_[pair_idx].end, split_arg.first + lagrange.order() + 1);

        lagrange.evaluate_derivative_table_at_x(split_arg.second);

        // == Expansion quantities ==

        const double innerprod =
            dot1.dipole().dot(e1.d0 * e0.d0 * dot0.dipole());
        const std::array<cmplx, 2> dyad{
            normalization(dr) * std::pow(c0, 2) *
                dot0.dipole().dot(e0.del_sq * e1.d0 * dot1.dipole()),
            normalization(dr) * std::pow(c0, 2) *
                dot1.dipole().dot(e1.del_sq * e0.d0 * dot0.dipole()),
        };

        for(int poly = 0; poly < lagrange.order() + 1; ++poly) {
          const int convolution_idx = split_arg.first + poly;
          const cmplx time =
              (lagrange.evaluations[2][poly] +
               2.0 * iu * omega_ * lagrange.evaluations[1][poly] -
               std::pow(omega_, 2) * lagrange.evaluations[0][poly]) *
              innerprod * normalization(dr);

          coefficients[pair_idx][convolution_idx][0] +=
              -time + dyad[0] * lagrange.evaluations[0][poly];

          if(pair.first == pair.second) continue;

          coefficients[pair_idx][convolution_idx][1] +=
              -time + dyad[1] * lagrange.evaluations[0][poly];
        }
      }
    }
  }

  return coefficients;
}
