#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
std::map<int, Detector> read_all_DRF_files (const std::string &DRF_prefix);

void Hires::iterate (std::vector<Sample> &samples)
{
  std::map<int, Detector> detectors (read_all_DRF_files (drf_prefix));

  Footprint footprints (radians_per_pix, ni, nj, min_sample_flux,
                        angle_tolerance, footprints_per_pix, detectors,
                        samples);
  wgt_image = footprints.calc_wgt_image (ni, nj);
  //
  // Standard hires/minimap processing - always happens.
  //
  arma::mat correction, correction_squared;
  footprints.compute_correction (ni, nj, flux_images[iteration - 1], iteration, 1,
                                 boost_func, boost_max_iter, correction,
                                 correction_squared);
  correction /= wgt_image;
  flux_images[iteration] = flux_images[iteration - 1] % correction; // Schur product

  arma::mat corr_sq_image = (correction_squared / wgt_image)
    - square (correction);
  cfv_images[iteration] = corr_sq_image;

  //
  // Optional HIRES beam generation
  //
  if (find (outfile_types.begin (), outfile_types.end (), "beam")
      != outfile_types.end ())
    {
      footprints.set_fluxes_to_sim_values (spike_image ());

      arma::mat correction, correction_squared;
      footprints.compute_correction (ni, nj, beam_images[iteration - 1], iteration,
                                     false, boost_func, boost_max_iter,
                                     correction, correction_squared);
      beam_images[iteration] = beam_images[iteration - 1] % correction
        / wgt_image; // Schur product
    }
  ++iteration;
}
}
