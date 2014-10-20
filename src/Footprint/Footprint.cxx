#include "../Footprint.hxx"

namespace hires
{
Footprint::Footprint (const double &radians_per_pix,
                      const std::array<int,2> &nxy,
                      const double &angle_tolerance,
                      const double &footprints_per_pix,
                      const std::map<int, Detector> &detectors,
                      const std::vector<Sample> &samples)
{
  size_t num_good = count_good_samples (radians_per_pix, nxy, samples);
  j0_im.reserve (num_good);
  j1_im.reserve (num_good);
  i0_im.reserve (num_good);
  i1_im.reserve (num_good);
  j0_ft.reserve (num_good);
  j1_ft.reserve (num_good);
  i0_ft.reserve (num_good);
  i1_ft.reserve (num_good);

  std::array<double,2> offset{{nxy[0] / 2.0, nxy[1] / 2.0}};

  for (size_t s = 0; s < samples.size (); ++s)
    {
      const size_t n (good[s].size ());
      for (size_t j = 0; j < n; ++j)
        {
          if (!good[s][j])
            continue;
          double xi ((samples[s].x[j] / radians_per_pix) + offset[0]),
              yi ((samples[s].y[j] / radians_per_pix) + offset[1]);

          int i_int (xi), j_int (yi);
          double i_frac (xi - i_int), j_frac (yi - j_int);

          responses.push_back (
              get_response (samples[s].id, i_frac, j_frac, samples[s].angle[j],
                            angle_tolerance, footprints_per_pix,
                            radians_per_pix, detectors));
          std::vector<int> bounds (
              compute_bounds (*(*responses.rbegin ()), i_int, j_int, nxy));

          signal.push_back (samples[s].signal[j]);
          j0_im.push_back (bounds[0]);
          j1_im.push_back (bounds[1]);
          i0_im.push_back (bounds[2]);
          i1_im.push_back (bounds[3]);
          j0_ft.push_back (bounds[4]);
          j1_ft.push_back (bounds[5]);
          i0_ft.push_back (bounds[6]);
          i1_ft.push_back (bounds[7]);
        }
    }
}
}
