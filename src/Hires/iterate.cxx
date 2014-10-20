#include <algorithm>

#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
void Hires::iterate (const bool &boosting)
{
  arma::mat correction;
  footprints.compute_correction (nxy, signal_image, iteration+1,
                                 boosting, boost_function, correction);
  // FIXME: Do it all in one loop
  correction /= weight_image;
  signal_image%=correction; // Schur product
  ++iteration;
}
}
