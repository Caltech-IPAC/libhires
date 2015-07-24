#include <fstream>
#include <boost/math/constants/constants.hpp>
#include <mlpack/methods/lars/lars.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "../Hires.hxx"

void hires::Hires::compute_elastic_net (const double &sigma_drf)
{
  Binned_Data binned_data(bin_data());
  arma::mat A(compute_response_function(sigma_drf,binned_data));
  
  const double data_scale=std::max(arma::max(arma::max(binned_data.data)),
                                   std::abs(arma::min(arma::min(binned_data.data))));
  std::cout << "variance: " << binned_data.variance << " "
            << std::sqrt(binned_data.variance) << " "
            << data_scale << "\n";

  mlpack::regression::LARS lars(true,binned_data.variance/data_scale,
                                binned_data.variance/(data_scale*data_scale));
  arma::vec lars_image;
  lars.Regress(A,binned_data.data,lars_image,false);

  elastic_net.set_size(nxy[0],nxy[1]);
  for (size_t ix=0; ix<nxy[0]; ++ix)
    for (size_t iy=0; iy<nxy[1]; ++iy)
      {
        elastic_net(iy,ix)=lars_image(ix+nxy[0]*iy);
      }
}

