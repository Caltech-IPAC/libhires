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
  std::map<std::tuple<int,int,int,int>,std::vector<double> > responses;

  std::vector<double> flux;
  std::vector<int> j0_im, j1_im, i0_im, i1_im, j0_ft, j1_ft, i0_ft, i1_ft;

  Footprint(const double &radians_per_pix, const int &NPIXi, const int &NPIXj,
            const double &min_sample_flux, const double &angle_tolerance,
            const double &footprints_per_pix,
            const std::map<int,Detector> &detectors,
            std::vector<Sample> &samples);

  double count_good_samples(const double &radians_per_pix,
                            const int &NPIXi, const int &NPIXj,
                            const std::vector<Sample> &samples);
  std::vector<double>
  get_response(const int &detector_id, const double &i_frac,
               const double &j_frac, const double &angle,
               const double &angle_tolerance,
               const double &footprints_per_pix,
               const std::map<int,Detector> &detectors);

  std::vector<double>
  generate_response(const int &detector_id, const double &i_offset,
                    const double &j_offset, const double &recomposed_angle,
                    const std::map<int,Detector> &detectors);
};

#endif
