#ifndef PULSE_INTERACTION_H
#define PULSE_INTERACTION_H

#include "../pulse.h"
#include "interaction.h"

class PulseInteraction : public InteractionBase {
 public:
  PulseInteraction(const std::shared_ptr<const DotVector>,
                   const std::shared_ptr<const Pulse>,
                   const double,
                   const double);
  virtual const ResultArray &evaluate(const int);
  virtual const ResultArray &evaluate_present_field(const int n)
  {
    return evaluate(n);
  }

 private:
  std::shared_ptr<const Pulse> pulse;
  const double hbar;
};

#endif
