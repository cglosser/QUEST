#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../src/integrator/history.h"
#include "../src/interactions/AIM/aim_interaction.h"
#include "../src/interactions/direct_interaction.h"
#include "../src/interactions/green_function.h"
#include "../src/math_utils.h"
#include "../src/quantum_dot.h"

BOOST_AUTO_TEST_SUITE(PROPAGATION)

struct Parameters {
  const double c0;
  const double dt;
  const int num_timesteps;
  const int interpolation_order;
  const int num_pts;

  const double max_time;

  Parameters()
      : c0{1},
        dt{1},
        num_timesteps{2048},
        interpolation_order{5},
        num_pts{2},
        max_time{num_timesteps * dt} {};
};

BOOST_AUTO_TEST_SUITE(RETARDED_GAUSSIAN)

struct SimulationStructures : Parameters {
  std::function<double(double)> src_fn, obs_fn;

  std::shared_ptr<DotVector> qds;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  SimulationStructures()
      : qds{std::make_shared<DotVector>()},
        history{std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_pts, num_timesteps / 2, num_timesteps)}
  {
    // == Initialize dots ==
    qds->emplace_back(Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Vector3d(0, 0, 1));
    qds->emplace_back(Eigen::Vector3d(10.5, 0.6, 0.7),
                      Eigen::Vector3d(0, 0, 1));

    const double dr = (qds->at(1).position() - qds->at(0).position()).norm();
    const double delay = dr / c0;

    const double mu = max_time / 2.0, sigma = max_time / 12.0;

    src_fn = [=](const double t) { return gaussian((t - mu) / sigma); };
    obs_fn = [=](const double t) { return src_fn(t - delay); };

    // == Initialize history ==
    for(int i = -num_timesteps / 2; i < num_timesteps; ++i) {
      history->array_[0][i][0] = Eigen::Vector2cd(0, src_fn(i * dt));
      history->array_[1][i][0] = Eigen::Vector2cd(0, 0);
    }
  }
};

BOOST_FIXTURE_TEST_CASE(DIRECT_INTERACTION, SimulationStructures)
{
  // == Initialize interaction ==
  Propagation::Identity<cmplx> identity_kernel;
  DirectInteraction direct_interaction(qds, history, identity_kernel,
                                       interpolation_order, c0, dt);

  // == Run propagation simulation ==
  for(int i = 0; i < num_timesteps; ++i) {
    BOOST_TEST_MESSAGE(i);

    const double obs_val_calculated = direct_interaction.evaluate(i)(1).real();
    const double obs_val_actual = obs_fn(i * dt);
    BOOST_REQUIRE_CLOSE_FRACTION(obs_val_calculated, obs_val_actual, 1e-9);
  }
}

BOOST_FIXTURE_TEST_CASE(AIM_INTERACTION, SimulationStructures)
{
  constexpr int expansion_order = 3;
  // == Initialize AIM structures ==
  Eigen::Vector3d spacing = Eigen::Vector3d::Ones() * c0 * dt;
  AIM::Grid grid(spacing, expansion_order, *qds);
  const int transit_steps =
      grid.max_transit_steps(c0, dt) + interpolation_order;

  Propagation::Identity<cmplx> identity_kernel;
  AIM::Interaction aim_interaction(
      qds, history, identity_kernel, spacing, interpolation_order,
      expansion_order, 1, c0, dt, AIM::Expansions::Retardation(transit_steps),
      AIM::Normalization::unit);

  // == Run propagation simulation ==
  for(int i = 0; i < num_timesteps; ++i) {
    BOOST_TEST_MESSAGE(i);

    const double obs_val_calculated = aim_interaction.evaluate(i)(1).real();
    const double obs_val_actual = obs_fn(i * dt);

    if(i > 16) {
      // Because of how AIM works the first several timesteps are required to
      // "warm up" the simulation and will likely have a large amount of error
      // but a negligible signal, thus this test simply ignores the first 16
      // steps. The tolerance is relatively large due to the identity kernel;
      // obs values don't scale as 1/r, so the grid approximation is much worse.
      // Of course, the approximation can be made much better with higher order
      // interpolations and/or expansions.
      BOOST_REQUIRE_CLOSE_FRACTION(obs_val_calculated, obs_val_actual, 1e-4);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()  // RETARDED_GAUSSIAN

BOOST_AUTO_TEST_SUITE_END()  // PROPAGATION
