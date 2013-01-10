#include <armadillo>

arma::mat create_spike_image(const int &n, const double &height,
                             const int &ni, const int &nj)
{
  arma::mat spike(nj,ni);
  spike.fill(0.000001);
  const int i_delta(ni/n), j_delta(nj/n);
  for(int i=i_delta/2;i<ni;i+=i_delta)
    for(int j=j_delta/2;j<nj;j+=j_delta)
      spike(j,i)=height;
  return spike;
}
