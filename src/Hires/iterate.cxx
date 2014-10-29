#include <algorithm>

#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
void Hires::iterate ()
{
  Eigen::MatrixXd correction(signal_image.rows(), signal_image.cols());
  footprints.compute_correction (signal_image, correction);
  signal_image+=correction;
  ++iteration;
}
}
