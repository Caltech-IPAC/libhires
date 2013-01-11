#include "../Footprint.hxx"

namespace hires
{
  double Footprint::count_good_samples(const double &radians_per_pix,
                                       const int &ni, const int &nj,
                                       const std::vector<Sample> &samples)
  {
    const double x_radius = (radians_per_pix * (ni-1)) / 2.0;
    const double y_radius = (radians_per_pix * (nj-1)) / 2.0;
    int total_samps = 0;
    int total_good  = 0;
    good.resize(samples.size());

    for(size_t i=0; i<samples.size(); ++i)
      {
        good[i].resize(samples[i].x.size());
        total_samps+=good[i].size();
        for(size_t j=0;j<good[i].size();++j)
          {
            good[i][j]=(samples[i].x[j]>-x_radius) && (samples[i].x[j]<x_radius)
              && (samples[i].y[j]>-y_radius) && (samples[i].y[j]<y_radius);
            if(good[i][j])
              ++total_good;
          }
      }

    LOG4CXX_INFO(logger, "image size degrees  = " << x_radius*2.0 << " x "
                 << y_radius*2.0 << "\n");

    LOG4CXX_INFO(logger, "IN data samples scanning complete\n");
    LOG4CXX_INFO(logger, total_samps << " data samples read\n");
    if(total_good==0)
      throw Exception("All data samples rejected (probably out of image)");
    LOG4CXX_INFO(logger, total_good << " data samples good\n");
    LOG4CXX_INFO(logger, total_samps-total_good << " data samples rejected\n");
    if(static_cast<double>(total_good)/total_samps < 0.5)
      LOG4CXX_WARN(logger, "More than 50% of data samples rejected\n");
    return total_good;
  }
}
