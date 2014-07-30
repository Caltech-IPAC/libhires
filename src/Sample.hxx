#pragma once

#include <valarray>
#include <utility>

namespace hires
{
class Sample
{
public:
  std::valarray<double> x, y, flux, angle;
  int id;

  Sample (const std::valarray<double> &X, const std::valarray<double> &Y,
          const std::valarray<double> &Flux, const int &Id,
          const std::valarray<double> &Angle)
      : x (X), y (Y), flux (Flux), angle (Angle), id (Id)
  {
  }
  Sample (const std::valarray<double> &X, const std::valarray<double> &Y,
          const std::valarray<double> &Flux, const int &Id)
      : Sample (X, Y, Flux, Id, std::valarray<double>())
  {
  }
  Sample () {}
};
}

