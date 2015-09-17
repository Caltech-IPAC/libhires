#include <mlpack/methods/lars/lars.hpp>
#include "../Hires.hxx"

void hires::Hires::compute_elastic_net (const double &sigma_drf)
{
  const size_t max_bins(64);
  Binned_Data binned_data (samples,nxy,radians_per_pix,max_bins);
  arma::mat A(compute_response_function(sigma_drf,binned_data));
  
  const double data_scale=std::max(arma::max(arma::max(binned_data.data)),
                                   std::abs(arma::min(arma::min(binned_data.data))));
  mlpack::regression::LARS lars(true,(binned_data.variance/data_scale)
                                *binned_data.num_bins*binned_data.num_bins/(nxy[0]*nxy[0]),
                                (binned_data.variance/(data_scale*data_scale))
                                *binned_data.num_bins*binned_data.num_bins/(nxy[0]*nxy[0]));
  arma::vec lars_image;
  lars.Regress(A,binned_data.data,lars_image,false);

  elastic_net.set_size(nxy[0],nxy[1]);
  for (size_t ix=0; ix<nxy[0]; ++ix)
    for (size_t iy=0; iy<nxy[1]; ++iy)
      {
        elastic_net(iy,ix)=lars_image(ix+nxy[0]*iy);
      }
}

