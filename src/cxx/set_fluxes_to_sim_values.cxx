/* Calculate simulated fluxes; replace fp.flux value
   Estimates what flux values would generate the input image */

#include "Footprint.hxx"
#include <armadillo>

void set_fluxes_to_sim_values(Footprint &fp,
                              const arma::mat &sim_image)
{
  for(size_t n=0;n<fp.flux.size();++n)
    {
      arma::mat integration(fp.j1_im[n]-fp.j0_im[n],fp.i1_im[n]-fp.i0_im[n]);
      for(int i=0;i<fp.i1_im[n]-fp.i0_im[n];++i)
        for(int j=0;j<fp.j1_im[n]-fp.j0_im[n];++j)
          {
            integration(j,i)=sim_image(j+fp.j0_im[n],i+fp.i0_im[n])
              *(*fp.responses[n])(j+fp.j0_ft[n],i+fp.i0_ft[n]);
          }
      fp.flux[n]=accu(integration);
    }
  LOG4CXX_INFO(logger,"Fluxes have been reset to simulated values\n");
}
