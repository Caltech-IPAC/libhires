#include "../Hires.hxx"

namespace hires
{
void Hires::compute_minimap ()
{
  double i_offset (nxy[0] / 2.0), j_offset (nxy[1] / 2.0);

  minimap.setZero ();
  hitmap.setZero ();

  for (size_t i = 0; i < samples.size (); ++i)
    {
      for (size_t j = 0; j < samples[i].x.size(); ++j)
        {
          double xi ((samples[i].x[j] / radians_per_pix) + i_offset);
          double yi ((samples[i].y[j] / radians_per_pix) + j_offset);
          int i_int (xi), j_int (yi);
          if(i_int<0 || i_int>=nxy[1] || j_int<0 || j_int>=nxy[0])
             continue;
          minimap (j_int, i_int) += samples[i].signal[j];
          hitmap (j_int, i_int) += 1.0;
        }
    }

  minimap = minimap.cwiseQuotient(hit_count);
}
}
