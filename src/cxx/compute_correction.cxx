#include <armadillo>
#include "Footprint.hxx"

void compute_correction(const int &nx, const int &ny,
                        const Footprint &fp,
                        const arma::mat &flux_image,
                        const int &iter, const bool &do_cfv,
                        const std::function<double (double)> &boost_func,
                        const int &boost_max_iter,
                        arma::mat &correction,
                        arma::mat &correction_squared)
{
  correction.zeros(nx,ny);
  if(do_cfv)
    correction_squared.zeros(nx,ny);

  std::string boost_mode;

  for(size_t n=0;n<fp.flux.size();++n)
    {
      arma::mat integration(fp.j1_im[n]-fp.j0_im[n],
                                  fp.i1_im[n]-fp.i0_im[n]);
      for(int i=0;i<fp.i1_im[n]-fp.i0_im[n];++i)
        for(int j=0;j<fp.j1_im[n]-fp.j0_im[n];++j)
          {
            integration(j,i)=flux_image(j+fp.j0_im[n],i+fp.i0_im[n])
              *(*fp.responses[n])(j+fp.j0_ft[n],i+fp.i0_ft[n]);
          }
      double flux_prime(accu(integration));
      double scale(fp.flux[n]/flux_prime);
      if(iter==1 && iter<=boost_max_iter)
        {
          scale=boost_func(scale);
          boost_mode="   (BOOSTED correction)\n";
        }
      for(int i=0;i<fp.i1_im[n]-fp.i0_im[n];++i)
        for(int j=0;j<fp.j1_im[n]-fp.j0_im[n];++j)
          {
            correction(j+fp.j0_im[n],i+fp.i0_im[n])+=
              (*fp.responses[n])(j+fp.j0_ft[n],i+fp.i0_ft[n])*scale;
          }
      if(do_cfv)
        {
          for(int i=0;i<fp.i1_im[n]-fp.i0_im[n];++i)
            for(int j=0;j<fp.j1_im[n]-fp.j0_im[n];++j)
              {
                correction_squared(j+fp.j0_im[n],i+fp.i0_im[n])+=
                  (*fp.responses[n])(j+fp.j0_ft[n],i+fp.i0_ft[n])*scale*scale;
              }
        }
    }
  LOG4CXX_INFO(logger,"Correction array compute, iter " << iter << "\n");
  if(!boost_mode.empty())
    LOG4CXX_INFO(logger,boost_mode);
}
