#include "pulse.h"

Pulse::Pulse(const double amplitude,
             const double delay,
             const double width,
             const double freq,
             const Eigen::Vector3d &wavevector,
             const Eigen::Vector3d &polarization,
             bool rotating)
    : amplitude(amplitude),
      delay(delay),
      width(width),
      freq(freq),
      wavevector(wavevector),
      polarization(polarization.normalized()),
      rotating(rotating)
{
}

Eigen::Vector3d Pulse::operator()(const Eigen::Vector3d &r,
                                  const double t) const
{
  const double arg = wavevector.dot(r) - freq * (t - delay);

  Eigen::Vector3d pulse_vector;
  pulse_vector = (amplitude / 2 * polarization) * gaussian(arg / width);

  pulse_vector *= rotating ? 1.0 : cos(arg);

  return pulse_vector;
}

std::ostream &operator<<(std::ostream &os, const Pulse &p)
{
  os << p.amplitude << " " << p.delay << " " << p.width << " " << p.freq << " "
     << p.wavevector.transpose() << " " << p.polarization.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, Pulse &p)
{
  is >> p.amplitude >> p.delay >> p.width >> p.freq >> p.wavevector[0] >>
      p.wavevector[1] >> p.wavevector[2] >> p.polarization[0] >>
      p.polarization[1] >> p.polarization[2];
  return is;
}

void set_reference_frame(
    Pulse &p, const bool rotating)  // there has to be a better way to do this
{
  p.rotating = rotating;
}

Pulse read_pulse_config(const std::string &filename, const bool rotating)
{
  std::ifstream ifs(filename);
  if(!ifs) throw std::runtime_error("Could not open " + filename);

  Pulse p;
  ifs >> p;
  set_reference_frame(p, rotating);

  return p;
}
