#include <algorithm>
#include <armadillo>
#include "../Hires.hxx"

namespace hires
{
arma::mat Hires::spike_image ()
{
  arma::mat spike (nxy[1], nxy[0]);
  spike.fill (0.000001);
  const int i_delta (std::max(nxy[0] / beam_spike_n,1)),
    j_delta (std::max(nxy[1] / beam_spike_n,1));
  for (int i = i_delta / 2; i < nxy[0]; i += i_delta)
    for (int j = j_delta / 2; j < nxy[1]; j += j_delta)
      spike (j, i) = beam_spike_height;
  return spike;
}
}
