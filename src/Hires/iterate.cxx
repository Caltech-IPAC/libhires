#include <algorithm>

#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
void Hires::iterate (const bool &boosting)
{
  arma::mat correction, correction_squared;
  footprints.compute_correction (nxy, signal_image, iteration+1, 1,
                                 boosting, boost_function, correction,
                                 correction_squared);
  // FIXME: Do it all in one loop
  correction /= weight_image;
  signal_image%=correction; // Schur product
  ++iteration;
}
}
