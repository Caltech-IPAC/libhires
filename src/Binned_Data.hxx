#pragma once

#include <array>

#include <armadillo>
#include "Sample.hxx"

namespace hires
{
class Binned_Data
{
public:
  arma::vec data;
  size_t num_bins;
  double variance;

  Binned_Data (const std::vector<Sample> &samples,
               const std::array<size_t,2> &nxy,
               const double &radians_per_pix,
               const size_t &max_bins);
};
}
