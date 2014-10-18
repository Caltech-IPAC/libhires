#include "../Footprint.hxx"

namespace hires
{
double Footprint::count_good_samples (const double &radians_per_pix,
                                      const std::array<int,2> &nxy,
                                      const std::vector<Sample> &samples)
{
  const std::array<double,2> radius{{(radians_per_pix * (nxy[0] - 1)) / 2.0,
        (radians_per_pix * (nxy[1] - 1)) / 2.0}};
  int total_samps = 0;
  int total_good = 0;
  good.resize (samples.size ());

  for (size_t i = 0; i < samples.size (); ++i)
    {
      good[i].resize (samples[i].x.size ());
      total_samps += good[i].size ();
      for (size_t j = 0; j < good[i].size (); ++j)
        {
          good[i][j] = (samples[i].x[j] > -radius[0])
                       && (samples[i].x[j] < radius[0])
                       && (samples[i].y[j] > -radius[1])
                       && (samples[i].y[j] < radius[1]);
          if (good[i][j])
            ++total_good;
        }
    }

  if (total_good == 0)
    throw Exception ("All data samples rejected (probably out of image)");
  if (static_cast<double>(total_good) / total_samps < 0.5)
    std::cerr << "Warning: Only " << total_good << " of " << total_samps
              << " data samples accepted\n";
  return total_good;
}
}
