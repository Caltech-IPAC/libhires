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

/// MCM with binning
void hires::Hires::compute_hires (const double &sigma_drf, const size_t &num_iterations)
{
  Binned_Data binned_data(bin_data());
  arma::mat A(compute_response_function(sigma_drf,binned_data));
  
  arma::vec f(nxy[0]*nxy[1]), f_new;
  const double total_flux=arma::mean(binned_data.data);
  f.fill(total_flux);
  for (size_t iteration=0; iteration<num_iterations; ++iteration)
    {
      arma::vec Af=A*f;
      arma::vec gAf=binned_data.data/Af;
      f_new=f%(A.t()*gAf);
      f=f_new*(total_flux/arma::mean(f_new));
      {
        std::ofstream outfile("mcm_" + std::to_string(iteration+1));
        for (size_t x=0; x<nxy[0]; ++x)
          for (size_t y=0; y<nxy[1]; ++y)
            outfile << x << " "
                    << y << " "
                    << f(x+nxy[0]*y) << " "
                    << "\n";
      }
    }
}

