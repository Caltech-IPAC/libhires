#include <algorithm>

#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
void Hires::iterate (const bool &boosting)
{
  Footprint footprints (radians_per_pix, nxy,
                        angle_tolerance, footprints_per_pix, detectors,
                        samples);
  wgt_image = footprints.calc_wgt_image (nxy);

  arma::mat correction, correction_squared;
  footprints.compute_correction (nxy, flux_images, iteration+1, 1,
                                 boosting, boost_function, correction,
                                 correction_squared);
  correction /= wgt_image;
  flux_images%=correction; // Schur product

  cfv_images = (correction_squared / wgt_image) - square (correction);

  if (output_types.find(Image_Type::hires_beam)!=output_types.end())
    {
      footprints.set_signals_to_sim_values (spike_image ());

      arma::mat correction, correction_squared;
      footprints.compute_correction (nxy, beam_images, iteration+1,
                                     false, boosting, boost_function,
                                     correction, correction_squared);
      beam_images%=correction;
      beam_images/=wgt_image; // Schur product
    }
  ++iteration;
}
}
