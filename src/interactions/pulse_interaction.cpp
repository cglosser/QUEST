#include "pulse_interaction.h"

PulseInteraction::PulseInteraction(const std::shared_ptr<const DotVector> &dots,
                                   const std::shared_ptr<const Pulse> pulse,
                                   const double hbar,
                                   const double dt)
    : Interaction(dots), pulse(std::move(pulse)), hbar(hbar), dt(dt)
{
}

const Interaction::ResultArray &PulseInteraction::evaluate(const int time_idx)
{
  const double time = time_idx * dt;

  for(size_t i = 0; i < dots->size(); ++i) {
    results[i] = (*pulse)((*dots)[i].position(), time);
  }

  return results;
}
