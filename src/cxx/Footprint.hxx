#ifndef HIRES_FOOTPRINT_HXX
#define HIRES_FOOTPRINT_HXX

#include <vector>
#include <valarray>
#include <map>
#include "Detector.hxx"
#include "Sample.hxx"

class Footprint
{
public:
  std::vector<std::valarray<bool> > good;

  std::vector<double> flux;
  std::vector<int> j0_im, j1_im, i0_im, i1_im, j0_ft, j1_ft, i0_ft, i1_ft;

  Footprint(const double &radians_per_pix, const int &NPIXi, const int &NPIXj,
            const double &min_sample_flux,
            const std::map<int,Detector> &detectors,
            std::vector<Sample> &samples);

  double count_good_samples(const double &radians_per_pix,
                            const int &NPIXi, const int &NPIXj,
                            const std::vector<Sample> &samples);
};

#endif
