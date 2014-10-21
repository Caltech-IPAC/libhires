#include <algorithm>

#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
void Hires::iterate (const bool &boosting)
{
  Eigen::MatrixXd correction;
  footprints.compute_correction (nxy, signal_image, iteration+1,
                                 boosting, boost_function, correction);
  signal_image=signal_image.cwiseProduct(correction.cwiseQuotient(weight_image)); // Schur product
  ++iteration;
}
}
