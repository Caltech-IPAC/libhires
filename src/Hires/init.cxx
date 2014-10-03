#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
std::map<int, Detector> read_all_DRF_files (const std::string &DRF_prefix);

void Hires::init (std::vector<Sample> &samples)
{
  std::map<int, Detector> detectors (read_all_DRF_files (drf_prefix));

  Footprint footprints (radians_per_pix, ni, nj, min_sample_flux,
                        angle_tolerance, footprints_per_pix, detectors,
                        samples);
  wgt_image = footprints.calc_wgt_image (ni, nj);
  // FIXME: Shouldn't I just use iteration?
  int iter_start;

  //
  // Standard hires/minimap processing - always happens.
  //
  flux_images = start_image (starting_image, iter_start);
  footprints.compute_minimap (radians_per_pix, ni, nj, samples, minimap,
                              hitmap);
  //
  // Optional HIRES beam generation
  //
  if (find (outfile_types.begin (), outfile_types.end (), "beam")
      != outfile_types.end ())
    {
      footprints.set_fluxes_to_sim_values (spike_image ());
      beam_images = start_image (beam_starting_image, iter_start);
    }
}
}
