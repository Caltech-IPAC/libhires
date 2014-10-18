/* Calculate simulated signals; replace fp.signal value
   Estimates what signal values would generate the input image */

#include "../Footprint.hxx"
#include <armadillo>

namespace hires
{
void Footprint::set_signals_to_sim_values (const arma::mat &sim_image)
{
  for (size_t n = 0; n < signal.size (); ++n)
    {
      arma::mat integration (j1_im[n] - j0_im[n], i1_im[n] - i0_im[n]);
      for (int i = 0; i < i1_im[n] - i0_im[n]; ++i)
        for (int j = 0; j < j1_im[n] - j0_im[n]; ++j)
          {
            integration (j, i) = sim_image (j + j0_im[n], i + i0_im[n])
                                 * (*responses[n])(j + j0_ft[n], i + i0_ft[n]);
          }
      signal[n] = accu (integration);
    }
}
}
