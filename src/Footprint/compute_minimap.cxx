#include <armadillo>
#include "../Footprint.hxx"
#include "../Sample.hxx"

namespace hires
{
void Footprint::compute_minimap (const double &radians_per_pix,
                                 const std::array<int,2> &nxy,
                                 const std::vector<Sample> &samples,
                                 arma::mat &minimap_image,
                                 arma::mat &minimap_hitmap) const
{
  double i_offset (nxy[0] / 2.0), j_offset (nxy[1] / 2.0);
  arma::mat integration (nxy[1], nxy[0]);
  arma::mat hitcnt (nxy[1], nxy[0]);

  integration.zeros ();
  hitcnt.zeros ();

  for (size_t i = 0; i < samples.size (); ++i)
    {
      for (size_t j = 0; j < good[i].size (); ++j)
        {
          if (!good[i][j])
            continue;

          double xi ((samples[i].x[j] / radians_per_pix) + i_offset);
          double yi ((samples[i].y[j] / radians_per_pix) + j_offset);
          int i_int (xi), j_int (yi);
          integration (j_int, i_int) += samples[i].signal[j];
          hitcnt (j_int, i_int) += 1.0;
        }
    }

  minimap_image = integration / hitcnt;
  minimap_hitmap = hitcnt;
}
}
