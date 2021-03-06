#include <cmath>

#include "../Hires.hxx"

namespace hires
{
// FIXME: This is largely superceded by bin_data.
void Hires::compute_minimap ()
{
  /// Offsets are nxy/2, not (nxy-1)/2, since the floor function needs
  /// an offset of 0.5.
  double i_offset (nxy[0] / 2.0), j_offset (nxy[1] / 2.0);

  minimap.set_size(nxy[0],nxy[1]);
  minimap.zeros();

  arma::mat hitmap;
  hitmap.set_size (nxy[0],nxy[1]);
  hitmap.zeros ();

  for (auto &sample: samples)
    {
      double xi ((sample.x / radians_per_pix) + i_offset);
      double yi ((sample.y / radians_per_pix) + j_offset);
      if(xi<0 || xi>=nxy[1] || yi<0 || yi>=nxy[0])
        continue;
      int i_int (std::floor(xi)), j_int (std::floor(yi));
      minimap (j_int, i_int) += sample.signal;
      hitmap (j_int, i_int) += 1.0;
    }
  minimap = minimap/hitmap;
}
}
