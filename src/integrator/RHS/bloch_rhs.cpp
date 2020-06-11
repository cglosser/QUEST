#include "bloch_rhs.h"

Integrator::BlochRHS::BlochRHS(
    const double dt,
    const std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history,
    std::vector<std::shared_ptr<InteractionBase>> interactions,
    std::vector<BlochFunctionType> rhs_functions)
    : Integrator::RHS<Eigen::Vector2cd>(dt, history),
      num_solutions(history->array_.shape()[0]),
      interactions(std::move(interactions)),
      rhs_functions(std::move(rhs_functions))
{
}

void Integrator::BlochRHS::evaluate_present(const int step) const
{
  auto eval_and_sum = [step](
                          const InteractionBase::ResultArray &r,
                          const std::shared_ptr<InteractionBase> &interaction) {
    return r + interaction->evaluate_present_field(step);
  };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

  auto projected_efields = std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);

  for(int solution = 0; solution < num_solutions; ++solution) {
    history->array_[solution][step][1] = rhs_functions[solution](
        history->array_[solution][step][0], projected_efields[solution]);
  }
}

void Integrator::BlochRHS::evaluate(const int step) const
{
  auto eval_and_sum = [step](
                          const InteractionBase::ResultArray &r,
                          const std::shared_ptr<InteractionBase> &interaction) {
    return r + interaction->evaluate(step);
  };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

  auto projected_efields = std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);

  for(int solution = 0; solution < num_solutions; ++solution) {
    history->array_[solution][step][1] = rhs_functions[solution](
        history->array_[solution][step][0], projected_efields[solution]);
  }
}
