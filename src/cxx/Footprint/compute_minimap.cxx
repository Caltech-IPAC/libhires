#include <armadillo>
#include "../Footprint.hxx"
#include "../Sample.hxx"

// Test for NaN with bitwise logic
static bool isNaN(float x) {
    return ((int&)x & 0x7fffffff) >= 0x7f800001;
}

namespace hires
{
  void Footprint::compute_minimap
  ( 
   const double &radians_per_pix,
   const int &nx, const int &ny,
   const std::vector<Sample> &samples,
   arma::mat &minimap_image,
   arma::mat &minimap_hitmap
  )
  {
    double i_offset(nx/2.0), j_offset(ny/2.0);
    arma::mat integration(nx, ny);
    arma::mat hitcnt(nx, ny);
  
    integration.zeros();
    hitcnt.zeros();

    for (size_t i=0; i < samples.size(); ++i)
      {
     
	for (size_t j=0; j < good[i].size(); ++j)
          { 
            if (!good[i][j] || isNaN(samples[i].flux[j])) continue;

	    double xi((samples[i].x[j]/radians_per_pix) + i_offset);
            double yi((samples[i].y[j]/radians_per_pix) + j_offset);
            int i_int(xi), j_int(yi);
         
	    integration(j_int,i_int) += samples[i].flux[j];
            hitcnt(j_int,i_int) += 1.0;
          }  
      }

      minimap_image = integration / hitcnt ;
      minimap_hitmap = hitcnt;
  }
}
