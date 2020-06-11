#include "interactions/AIM/direct.h"

AIM::DirectInteraction::DirectInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    Propagation::Kernel<cmplx> &kernel,
    const int interp_order,
    const double c0,
    const double dt,
    std::shared_ptr<const std::vector<Grid::ipair_t>> interaction_pairs)
    : HistoryInteraction(dots, history, interp_order, c0, dt),
      interaction_pairs_{std::move(interaction_pairs)},
      shape_(
          {{static_cast<int>(interaction_pairs_->size()), interp_order + 1}}),
      floor_delays_(shape_[0]),
      coefficients_(coefficient_table(kernel, floor_delays_))
{
}

const InteractionBase::ResultArray &AIM::DirectInteraction::evaluate(
    const int time_idx)
{
  constexpr int RHO_01 = 1;

  results.setZero();
  past_terms_of_convolution.setZero();

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];

    for(int i = 1; i < shape_[1]; ++i) {
      const int s =
          std::max(time_idx - floor_delays_[pair_idx] - i,
                   static_cast<int>(history->array_.index_bases()[1]));

      past_terms_of_convolution[pair.first] +=
          (history->array_[pair.second][s][0])[RHO_01] *
          coefficients_[pair_idx][i];
      past_terms_of_convolution[pair.second] +=
          (history->array_[pair.first][s][0])[RHO_01] *
          coefficients_[pair_idx][i];
    }

    const int s = std::max(time_idx - floor_delays_[pair_idx],
                           static_cast<int>(history->array_.index_bases()[1]));

    results[pair.first] += (history->array_[pair.second][s][0])[RHO_01] *
                           coefficients_[pair_idx][0];
    results[pair.second] += (history->array_[pair.first][s][0])[RHO_01] *
                            coefficients_[pair_idx][0];
  }
  results += past_terms_of_convolution;
  return results;
}

const InteractionBase::ResultArray &
AIM::DirectInteraction::evaluate_present_field(const int time_idx)
{
  constexpr int RHO_01 = 1;

  results.setZero();

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];

    const int s = std::max(time_idx - floor_delays_[pair_idx],
                           static_cast<int>(history->array_.index_bases()[1]));

    results[pair.first] += (history->array_[pair.second][s][0])[RHO_01] *
                           coefficients_[pair_idx][0];
    results[pair.second] += (history->array_[pair.first][s][0])[RHO_01] *
                            coefficients_[pair_idx][0];
  }
  results += past_terms_of_convolution;
  return results;
}

boost::multi_array<cmplx, 2> AIM::DirectInteraction::coefficient_table(
    Propagation::Kernel<cmplx> &kernel, std::vector<int> &floor_delays) const
{
  boost::multi_array<cmplx, 2> coefs(shape_);
  std::fill(coefs.data(), coefs.data() + coefs.num_elements(), cmplx(0.0, 0.0));
  floor_delays.resize(shape_[0]);

  Interpolation::UniformLagrangeSet lagrange(shape_[1] - 1);

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];

    if(pair.first == pair.second) continue;

    Eigen::Vector3d dr(separation((*dots)[pair.first], (*dots)[pair.second]));
    auto delay = split_double(dr.norm() / (c0 * dt));

    floor_delays[pair_idx] = delay.first;
    lagrange.evaluate_derivative_table_at_x(delay.second, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        kernel.coefficients(dr, lagrange));

    for(int i = 0; i < shape_[1]; ++i) {
      coefs[pair_idx][i] = (*dots)[pair.second].dipole().dot(
          interp_dyads[i] * (*dots)[pair.first].dipole());
    }
  }

  return coefs;
}
