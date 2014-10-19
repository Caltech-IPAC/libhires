#include "Detector.hxx"

namespace hires
{
std::map<int, Detector> read_DRF (const boost::filesystem::path &DRF_file)
{
  std::map<int, Detector> detectors;
  Detector d (DRF_file);
  /* A trick so that the underlying array is not copied */
  std::swap (detectors[d.id], d);
  return detectors;
}
}
