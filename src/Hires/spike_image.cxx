#include <armadillo>
#include "../Hires.hxx"

namespace hires
{
arma::mat Hires::spike_image ()
{
  arma::mat spike (nj, ni);
  spike.fill (0.000001);
  const int i_delta (ni / beam_spike_n), j_delta (nj / beam_spike_n);
  for (int i = i_delta / 2; i < ni; i += i_delta)
    for (int j = j_delta / 2; j < nj; j += j_delta)
      spike (j, i) = beam_spike_height;
  return spike;
}
}
