#include <boost/math/constants/constants.hpp>
#include "../Hires.hxx"

arma::mat hires::Hires::compute_response_function
(const double &sigma_drf,
 const Binned_Data &binned_data)
{
  const size_t bins=binned_data.num_bins;
  const double sigma_drf2=sigma_drf*sigma_drf;
  const double min_x(0), min_y(0), delta_bin_x(nxy[0]*radians_per_pix/bins),
    delta_bin_y(nxy[1]*radians_per_pix/bins), pi(boost::math::constants::pi<double>());
  
  arma::mat response(bins*bins,nxy[0]*nxy[1]);
  for(size_t iy=0; iy<nxy[1]; ++iy)
    {
      const double y=min_y + (iy+0.5)*radians_per_pix;
      for(size_t ix=0; ix<nxy[0]; ++ix)
        {
          const double x=min_x + (ix+0.5)*radians_per_pix;

          for(size_t iy_bin=0; iy_bin<bins; ++iy_bin)
            {
              const double bin_y=min_y + (iy_bin+0.5)*delta_bin_y;
              for(size_t ix_bin=0; ix_bin<bins; ++ix_bin)
                {
                  const double bin_x=min_x + (ix_bin+0.5)*delta_bin_x;
                  double dx(x - bin_x), dy(y-bin_y);
                  response(ix_bin + bins*iy_bin,ix + nxy[0]*iy)=
                    radians_per_pix*radians_per_pix
                    *exp(-(dx*dx+dy*dy)/(2*sigma_drf2))/(2*pi*sigma_drf2);
                }
            }
        }
    }
  return response;
}

