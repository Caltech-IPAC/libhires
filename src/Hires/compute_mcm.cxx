#include "../Hires.hxx"

/// MCM with binning
void hires::Hires::compute_mcm (const double &sigma_drf,
                                const size_t &num_iterations)
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
    }
  hires.set_size(nxy[0],nxy[1]);
  for (size_t ix=0; ix<nxy[0]; ++ix)
    for (size_t iy=0; iy<nxy[1]; ++iy)
      {
        hires(iy,ix)=f(ix+nxy[0]*iy);
      }
}

