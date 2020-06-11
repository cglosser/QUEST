#ifndef INTERACTION_H
#define INTERACTION_H

#include <Eigen/Dense>
#include <memory>

#include "../quantum_dot.h"

class InteractionBase {
 public:
  typedef Eigen::Array<cmplx, Eigen::Dynamic, 1> ResultArray;

  InteractionBase(const std::shared_ptr<const DotVector> dots, const double dt)
      : dots(std::move(dots)),
        results(dots ? dots->size() : 0),
        past_terms_of_convolution(dots ? dots->size() : 0),
        dt(dt){};
  // const cmplx &operator[](const int i) const { return results[i]; }
  virtual const ResultArray &evaluate(const int) = 0;
  virtual const ResultArray &evaluate_present_field(const int) = 0;

  virtual ~InteractionBase(){};

 protected:
  std::shared_ptr<const DotVector> dots;
  ResultArray results;
  ResultArray past_terms_of_convolution;
  double dt;
};

#endif
