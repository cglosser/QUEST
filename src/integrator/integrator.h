#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <string>

#include "RHS/rhs.h"
#include "history.h"
#include "logging.h"
#include "math_utils.h"
#include "weights.h"

namespace Integrator {
  template <class soltype>
  class PredictorCorrector;
}

template <class soltype>
class Integrator::PredictorCorrector {
 public:
  PredictorCorrector(const double,
                     const int,
                     const int,
                     const int,
                     const double,
                     const std::shared_ptr<Integrator::History<soltype>>,
                     std::unique_ptr<Integrator::RHS<soltype>>);
  void solve(const log_level_t = log_level_t::LOG_NOTHING) const;

 private:
  int num_solutions, time_idx_ubound, num_corrector_steps;
  double dt;
  Weights weights;
  std::shared_ptr<Integrator::History<soltype>> history;
  std::unique_ptr<Integrator::RHS<soltype>> rhs;

  void solve_step(const int) const;
  void predictor(const int) const;
  void corrector(const int) const;

  void log_percentage_complete(const int) const;
};

template <class soltype>
Integrator::PredictorCorrector<soltype>::PredictorCorrector(
    const double dt,
    const int num_corrector_steps,
    const int n_lambda,
    const int n_time,
    const double radius,
    const std::shared_ptr<Integrator::History<soltype>> history,
    std::unique_ptr<Integrator::RHS<soltype>> rhs)
    : num_solutions(history->array_.shape()[0]),
      time_idx_ubound(history->array_.index_bases()[1] +
                      history->array_.shape()[1]),
      dt(dt),
      num_corrector_steps(num_corrector_steps),
      weights(n_lambda, n_time, radius),
      history(std::move(history)),
      rhs(std::move(rhs))
{
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::solve(
    const log_level_t log_level) const
{
  for(int step = 0; step < time_idx_ubound; ++step) {
    solve_step(step);
    if(log_level >= log_level_t::LOG_INFO) log_percentage_complete(step);
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::solve_step(const int step) const
{
  assert(0 <= step && step < time_idx_ubound);

  predictor(step);
  rhs->evaluate(step);

  for(int m = 0; m < num_corrector_steps; ++m) {
    corrector(step);
    rhs->evaluate_present(step);
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::predictor(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history->array_[sol_idx][step][0] +=
          history->array_[sol_idx][start + h][0] * weights.ps(0, h) +
          history->array_[sol_idx][start + h][1] * weights.ps(1, h) * dt;
    }
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::corrector(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    history->array_[sol_idx][step][0] =
        weights.future_coef * history->array_[sol_idx][step][1] * dt;
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history->array_[sol_idx][step][0] +=
          history->array_[sol_idx][start + h][0] * weights.cs(0, h) +
          history->array_[sol_idx][start + h][1] * weights.cs(1, h) * dt;
    }
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::log_percentage_complete(
    const int step) const
{
  if(step % (time_idx_ubound / 10) == 0) {
    std::cout << "\t" << static_cast<int>(10.0 * step / time_idx_ubound)
              << std::endl;
  }
}

#endif
