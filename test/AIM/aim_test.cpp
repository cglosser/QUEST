#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/AIM"

BOOST_AUTO_TEST_SUITE(AIM)

template <int num_dots>
struct PARAMETERS {
  int interpolation_order, expansion_order, num_steps;
  double c, dt, total_time;

  Eigen::Array3d spacing;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  PARAMETERS()
      : interpolation_order(4),
        expansion_order(1),
        num_steps(1024),

        c(1),
        dt(1),
        total_time(1024),

        spacing(Eigen::Array3d(1, 1, 1) * c * dt),

        dots(std::make_shared<DotVector>(num_dots)),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, 10, num_steps))
  {
    std::cout.precision(17);

    history->fill(Eigen::Vector2cd::Zero());
    for(int t = -10; t < num_steps; ++t) {
      history->array_[0][t][0](RHO_01) = src(t * dt);
    }
  };

  double src(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 6.0);
    return gaussian(arg);
  }
};

BOOST_AUTO_TEST_SUITE(IDENTITY_KERNEL)

BOOST_FIXTURE_TEST_CASE(DIRECT_1, PARAMETERS<2>)
{
  dots->at(0) = QuantumDot({0.5, 0.5, 0.5}, {0, 0, 1});
  dots->at(1) = QuantumDot({0.5, 0.5, 1.5}, {0, 0, 1});
  AIM::Grid grid(spacing, dots, expansion_order);

  Eigen::Vector3d dr = dots->at(1).position() - dots->at(0).position();

  Propagation::Identity<cmplx> gf;
  AIM::DirectInteraction direct1(dots, history, gf, interpolation_order, c, dt,
                                 grid);

  for(int t = 0; t < num_steps; ++t) {
    const auto prop = direct1.evaluate(t);
    BOOST_CHECK_EQUAL(std::norm(prop(0)), 0);
    BOOST_CHECK_SMALL(std::norm(prop(1) - src(t * dt - dr.norm() / c)), 1e-14);
  }
}

BOOST_FIXTURE_TEST_CASE(DIRECT_2, PARAMETERS<2>)
{
  dots->at(0) = QuantumDot({0.5, 0.5, 0.5}, {0, 0, 1});
  dots->at(1) = QuantumDot({0.5, 0.5, 2.5}, {0, 0, 1});
  AIM::Grid grid(spacing, dots, expansion_order);

  Eigen::Vector3d dr = dots->at(1).position() - dots->at(0).position();

  Propagation::Identity<cmplx> gf;
  AIM::DirectInteraction direct1(dots, history, gf, interpolation_order, c, dt,
                                 grid);

  for(int t = 0; t < num_steps; ++t) {
    const auto prop = direct1.evaluate(t);
    BOOST_CHECK_EQUAL(std::norm(prop(0)), 0);
    BOOST_CHECK_SMALL(std::norm(prop(1) - src(t * dt - dr.norm() / c)), 1e-14);
  }
}

BOOST_FIXTURE_TEST_CASE(DIRECT_NIL, PARAMETERS<2>)
{
  dots->at(0) = QuantumDot({0.5, 0.5, 0.5}, {0, 0, 1});
  dots->at(1) = QuantumDot({0.5, 0.5, 3.5}, {0, 0, 1});
  AIM::Grid grid(spacing, dots, expansion_order);

  Propagation::Identity<cmplx> gf;
  AIM::DirectInteraction direct1(dots, history, gf, interpolation_order, c, dt,
                                 grid);

  for(int t = 0; t < num_steps; ++t) {
    const auto prop = direct1.evaluate(t);
    BOOST_CHECK_EQUAL(std::norm(prop(0)), 0);
    BOOST_CHECK_SMALL(std::norm(prop(1)), 1e-14);
  }
}

BOOST_AUTO_TEST_SUITE_END()  // IDENTITY_KERNEL

BOOST_AUTO_TEST_SUITE_END()  // AIM
