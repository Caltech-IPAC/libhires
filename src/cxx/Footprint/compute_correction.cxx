#include <armadillo>
#include "../Footprint.hxx"

// Test for NaN with bitwise logic
bool isNaN(float x) {
    return ((int&)x & 0x7fffffff) >= 0x7f800001;
}

namespace hires
{
  void Footprint::compute_correction
  (const int &nx, const int &ny,
   const arma::mat &flux_image,
   const int &iter, const bool &do_cfv,
   const std::function<double (double)> &boost_func,
   const int &boost_max_iter,
   arma::mat &correction,
   arma::mat &correction_squared)
  {
    correction.zeros(ny,nx);
    if(do_cfv)
      correction_squared.zeros(ny,nx);

    for(size_t n=0;n<flux.size();++n)
       {
// Don't process TOI samples having NaN values.
        if (isNaN(flux[n])) continue;

        arma::mat integration(j1_im[n]-j0_im[n],
                              i1_im[n]-i0_im[n]);
        for(int i=0;i<i1_im[n]-i0_im[n];++i)
          for(int j=0;j<j1_im[n]-j0_im[n];++j)
            {
//               if (!isNaN(flux_image(j+j0_im[n],i+i0_im[n]))) {
                  integration(j,i)=flux_image(j+j0_im[n],i+i0_im[n])
                      *(*responses[n])(j+j0_ft[n],i+i0_ft[n]);
//               } 
            } 

        double flux_prime(accu(integration));

        double scale(flux[n]/flux_prime);
        if(iter==1 && iter<=boost_max_iter)
          {
            scale=boost_func(scale);
          }
        for(int i=0;i<i1_im[n]-i0_im[n];++i)
          for(int j=0;j<j1_im[n]-j0_im[n];++j)
            {
              correction(j+j0_im[n],i+i0_im[n])+=
                (*responses[n])(j+j0_ft[n],i+i0_ft[n])*scale;
            }

        if(do_cfv)
          {
            for(int i=0;i<i1_im[n]-i0_im[n];++i)
              for(int j=0;j<j1_im[n]-j0_im[n];++j)
                {
                  correction_squared(j+j0_im[n],i+i0_im[n])+=
                    (*responses[n])(j+j0_ft[n],i+i0_ft[n])*scale*scale;
                }
          }
      }
  }
}
