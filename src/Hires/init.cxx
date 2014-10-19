#include <algorithm>

#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
void Hires::init ()
{
  if(output_types.find(Image_Type::minimap_image)!=output_types.end()
     || output_types.find(Image_Type::minimap_hitmap)!=output_types.end())
    compute_minimap ();

  if(running_hires())
    {
      Footprint footprints (radians_per_pix, nxy, min_sample_flux,
                            angle_tolerance, footprints_per_pix, detectors,
                            samples);
      wgt_image = footprints.calc_wgt_image (nxy);
      // FIXME: Shouldn't I just use iteration?
      int iter_start;

      flux_images = start_image (starting_image, iter_start);

      if (output_types.find(Image_Type::hires_beam)!=output_types.end())
        {
          footprints.set_signals_to_sim_values (spike_image ());
          beam_images = start_image (beam_starting_image, iter_start);
        }
    }
}
}
