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

  {
    std::ofstream outfile("elastic_net");
    for (size_t x=0; x<nxy[0]; ++x)
      for (size_t y=0; y<nxy[1]; ++y)
        outfile << x << " "
                << y << " "
                << lars_image(x+nxy[0]*y) << " "
                << "\n";
  }
}

