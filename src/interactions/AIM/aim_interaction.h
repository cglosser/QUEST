#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include "interactions/AIM/direct.h"
#include "interactions/AIM/farfield.h"
#include "interactions/AIM/nearfield.h"

namespace AIM {
  class Interaction;
}

class AIM::Interaction final : public InteractionBase {
 public:
  Interaction(const std::shared_ptr<DotVector> dots,
              const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>
                  history,
              Propagation::Kernel<cmplx> &kernel,
              const Eigen::Array3d &spacing,
              const int interp_order,
              const int expansion_order,
              const int border,
              const double c0,
              const double dt,
              Expansions::ExpansionFunction expansion_function,
              Normalization::SpatialNorm normalization,
              const double omega = 0)
      : InteractionBase(dots, dt),
        grid{std::make_shared<Grid>(spacing, expansion_order, *dots)},
        expansion_table{std::make_shared<Expansions::ExpansionTable>(
            Expansions::LeastSquaresExpansionSolver::get_expansions(
                expansion_order, *grid, *dots))},
        nearfield_pairs{std::make_shared<std::vector<Grid::ipair_t>>(
            grid->nearfield_point_pairs(border, *dots))},

        ff(dots,
           history,
           interp_order,
           c0,
           dt,
           grid,
           expansion_table,
           expansion_function,
           normalization),

        nf(dots,
           history,
           5,
           c0,
           dt,
           grid,
           expansion_table,
           nullptr,
           normalization,
           nearfield_pairs,
           omega),

        direct(dots, history, kernel, interp_order, c0, dt, nearfield_pairs)
  {
  }

  const ResultArray &evaluate(const int t)
  {
    results =
        (ff.evaluate(t).conjugate() - nf.evaluate(t)) + direct.evaluate(t);
    return results;
  }

  const ResultArray &evaluate_present_field(const int t)
  {
    // I DON'T KNOW WHY THAT NEEDS A CONJUGATE!!!
    results = (ff.evaluate_present_field(t).conjugate() -
               nf.evaluate_present_field(t)) +
              direct.evaluate_present_field(t);
    return results;
  }

 private:
  std::shared_ptr<Grid> grid;
  std::shared_ptr<Expansions::ExpansionTable> expansion_table;
  std::shared_ptr<std::vector<Grid::ipair_t>> nearfield_pairs;

  Farfield ff;
  Nearfield nf;
  DirectInteraction direct;
};

#endif
