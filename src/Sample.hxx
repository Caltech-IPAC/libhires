#pragma once

#include <valarray>
#include <utility>

namespace hires
{
class Sample
{
public:
  double x, y, signal, angle;
  Sample(const double &X, const double &Y, const double &Signal,
         const double &Angle):
    x(X), y(Y), signal(Signal), angle(Angle) {}
};
}

