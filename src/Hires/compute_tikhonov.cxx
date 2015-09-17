#include <mlpack/methods/linear_regression/linear_regression.hpp>
#include "../Hires.hxx"

void hires::Hires::compute_tikhonov (const double &sigma_drf)
{
  // const size_t max_bins(2*nxy[0]*radians_per_pix/sigma_drf);
  const size_t max_bins(128);
  Binned_Data binned_data (samples,nxy,radians_per_pix,max_bins);
  arma::mat A(compute_response_function(sigma_drf,binned_data));
  
  const double data_scale=std::max(arma::max(arma::max(binned_data.data)),
                                   std::abs(arma::min(arma::min(binned_data.data))));
  std::cout << "scale: " << data_scale << " " << binned_data.variance << "\n";
  mlpack::regression::LinearRegression
    regression(A.t(),binned_data.data,(binned_data.variance/(data_scale*data_scale))
               *binned_data.num_bins*binned_data.num_bins/(nxy[0]*nxy[0]));
    // regression(A.t(),binned_data.data,binned_data.num_bins/sqrt(nxy[0]*nxy[1]));
    // regression(A.t(),binned_data.data,1);
    // regression(A.t(),binned_data.data,binned_data.variance/data_scale);
  arma::vec image=regression.Parameters ();

  tikhonov.set_size(nxy[0],nxy[1]);
  for (size_t ix=0; ix<nxy[0]; ++ix)
    for (size_t iy=0; iy<nxy[1]; ++iy)
      {
        /// The tikhonov image has an offset in the first spot.
        tikhonov(iy,ix)=image(1+ix+nxy[0]*iy)+image(0);
      }
}

