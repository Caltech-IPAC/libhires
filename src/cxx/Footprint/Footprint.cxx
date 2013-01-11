#include <log4cxx/logger.h>
#include <log4cxx/fileappender.h>

#include "../Footprint.hxx"

namespace hires
{
  Footprint::Footprint(const double &radians_per_pix, const int &ni,
                       const int &nj,
                       const double &min_sample_flux,
                       const double &angle_tolerance,
                       const double &footprints_per_pix,
                       const std::map<int,Detector> &detectors,
                       std::vector<Sample> &samples)
  {
    size_t num_good=count_good_samples(radians_per_pix,ni,nj,samples);
    j0_im.reserve(num_good);
    j1_im.reserve(num_good);
    i0_im.reserve(num_good);
    i1_im.reserve(num_good);
    j0_ft.reserve(num_good);
    j1_ft.reserve(num_good);
    i0_ft.reserve(num_good);
    i1_ft.reserve(num_good);
  
    double i_offset(ni/2.0), j_offset(nj/2.0);

    int n_fluxes_reset(0);
    int num_footprints(0);

    for(size_t i=0;i<samples.size();++i)
      {
        const size_t n(good[i].size());

        /* Compute angles if needed */
        if(samples[i].angle.size()==0)
          {
            samples[i].angle.resize(n);
            double max_delta=radians_per_pix*16;
            for(size_t j=0;j<n;++j)
              {
                if(!good[i][j])
                  continue;
                if(j==n-1)
                  {
                    samples[i].angle[j]=samples[i].angle[j-1];
                    good[i][j]=good[i][j-1];
                  }
                else
                  {
                    double x_delta(samples[i].x[j+1]-samples[i].x[j]),
                      y_delta(samples[i].y[j+1]-samples[i].y[j]);
                    if(std::abs(x_delta) + std::abs(y_delta)>max_delta)
                      {
                        /* non-continuous */
                        if(j>0)
                          {
                            samples[i].angle[j]=samples[i].angle[j-1];
                            good[i][j]=good[i][j-1];
                          }
                        else
                          {
                            good[i][j]=false; /* Discard sample if we
                                                 can not compute
                                                 angle */
                          }
                      }
                    else
                      {
                        samples[i].angle[j]=atan2(y_delta,x_delta);
                      }
                  }
              }
          }

        for(size_t j=0;j<n;++j)
          {
            if(!good[i][j])
              continue;
            double xi((samples[i].x[j]/radians_per_pix)+i_offset),
              yi((samples[i].y[j]/radians_per_pix)+j_offset);

            int i_int(xi), j_int(yi);
            double i_frac(xi-i_int), j_frac(yi-j_int);

            responses.push_back(get_response(samples[i].id,i_frac,j_frac,
                                             samples[i].angle[j],
                                             angle_tolerance,
                                             footprints_per_pix,
                                             radians_per_pix,
                                             detectors));
            std::vector<int> bounds(compute_bounds(*(*responses.rbegin()),
                                                   i_int,j_int,ni,nj));

            if(samples[i].flux[j]<min_sample_flux)
              {
                samples[i].flux[j]=min_sample_flux;
                ++n_fluxes_reset;
              }
            flux.push_back(samples[i].flux[j]);
            j0_im.push_back(bounds[0]);
            j1_im.push_back(bounds[1]);
            i0_im.push_back(bounds[2]);
            i1_im.push_back(bounds[3]);
            j0_ft.push_back(bounds[4]);
            j1_ft.push_back(bounds[5]);
            i0_ft.push_back(bounds[6]);
            i1_ft.push_back(bounds[7]);
            ++num_footprints;
          }
      }
    LOG4CXX_INFO(logger,"Footprint creation complete\n");
    LOG4CXX_INFO(logger,num_footprints << " footprints created\n");
    LOG4CXX_INFO(logger,"(" << responses_complete.size()
                 << " full response arrays created)\n");
    LOG4CXX_INFO(logger,n_fluxes_reset << " sample fluxes reset to minimum\n");
  }
}
