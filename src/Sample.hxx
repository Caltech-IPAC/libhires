#pragma once

#include <valarray>
#include <utility>

namespace hires
{
class Sample
{
public:
  std::valarray<double> x, y, signal, angle;
  int id;

  Sample (const std::valarray<double> &X, const std::valarray<double> &Y,
          const std::valarray<double> &Signal, const int &Id,
          const std::valarray<double> &Angle)
      : x (X), y (Y), signal (Signal), angle (Angle), id (Id)
  {
  }
  Sample (const std::valarray<double> &X, const std::valarray<double> &Y,
          const std::valarray<double> &Signal, const int &Id)
      : Sample (X, Y, Signal, Id, std::valarray<double>())
  {
  }
  Sample () {}
};
}

