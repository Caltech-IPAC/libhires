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

  Footprint(const double &radians_per_pix, const int &NPIXi, const int &NPIXj,
            const std::map<int,Detector> &detectors,
            const std::vector<Sample> &samples)
  {
    size_t total_good=count_good_samples(radians_per_pix,NPIXi,NPIXj,samples);
  }

  double count_good_samples(const double &radians_per_pix,
                            const int &NPIXi, const int &NPIXj,
                            const std::vector<Sample> &samples);
};

#endif
