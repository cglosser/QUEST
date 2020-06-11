#ifndef ODE_RHS_H
#define ODE_RHS_H

#include <cmath>
#include <functional>
#include <vector>
#include "rhs.h"

namespace Integrator {
  class ODE_RHS;
}

class Integrator::ODE_RHS : public RHS<double> {
 public:
  ODE_RHS(const double,
          const std::shared_ptr<Integrator::History<double>>,
          const std::vector<std::function<double(double, double)>>);
  void evaluate_present(const int) const;
  void evaluate(const int) const;

 private:
  std::vector<std::function<double(double, double)>> rhs_functions;
};

Integrator::ODE_RHS::ODE_RHS(
    const double dt,
    const std::shared_ptr<Integrator::History<double>> history,
    const std::vector<std::function<double(double, double)>> rhs_functions)
    : RHS(dt, history), rhs_functions(std::move(rhs_functions))
{
}

void Integrator::ODE_RHS::evaluate(const int n) const
{
  const double time = n * dt;
  for(int i = 0; i < static_cast<int>(history->array_.shape()[0]); ++i) {
    history->array_[i][n][1] = rhs_functions[i](history->array_[i][n][0], time);
  }
}

void Integrator::ODE_RHS::evaluate_present(const int n) const { evaluate(n); }

#endif
