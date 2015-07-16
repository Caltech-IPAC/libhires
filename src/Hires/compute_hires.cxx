#include "../Hires.hxx"
#include "Detector.hxx"
#include "Footprint.hxx"

std::map<int, Detector> hires::read_DRF (const boost::filesystem::path &DRF_file);


hires::Hires::compute_hires (const boost::filesystem::path &Drf_file)
{
  std::map<int, Detector> detectors (read_DRF (drf_file));
  // FIXME: What use is this??
  int footprints_per_pix (1);
  Footprint footprints (radians_per_pix, nxy, angle_tolerance,
                        footprints_per_pix, detectors, samples);
  hires_image.resize(nxy[0],nxy[1]);
  hires_image.setZero();
}

