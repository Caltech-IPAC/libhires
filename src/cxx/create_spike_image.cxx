#include <armadillo>

arma::mat create_spike_image(const int &n, const double &height,
                             const int &NPIXi, const int &NPIXj)
{
  arma::mat spike(NPIXj,NPIXi);
  spike.fill(0.000001);
  const int i_delta(NPIXi/n), j_delta(NPIXj/n);
  for(int i=i_delta/2;i<NPIXi;i+=i_delta)
    for(int j=j_delta/2;j<NPIXj;j+=j_delta)
      spike(j,i)=height;
  return spike;
}
