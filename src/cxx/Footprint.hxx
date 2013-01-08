#ifndef HIRES_FOOTPRINT_HXX
#define HIRES_FOOTPRINT_HXX

#include <vector>
#include <valarray>
#include <map>
#include <Eigen/Dense>
#include "Detector.hxx"
#include "Sample.hxx"

class Footprint
{
public:
  std::vector<std::valarray<bool> > good;
  std::map<std::tuple<int,int,int,int>,Eigen::MatrixXd> responses_complete;
  std::vector<Eigen::MatrixXd> responses;

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
  Eigen::MatrixXd
  get_response(const int &detector_id, const double &i_frac,
               const double &j_frac, const double &angle,
               const double &angle_tolerance,
               const double &footprints_per_pix,
               const std::map<int,Detector> &detectors);

  Eigen::MatrixXd
  generate_response(const int &detector_id, const double &i_offset,
                    const double &j_offset, const double &recomposed_angle,
                    const std::map<int,Detector> &detectors);

  std::vector<int> compute_bounds(const Eigen::MatrixXd &response,
                                  const int &i_center, const int &j_center,
                                  const int &NPIXi, const int &NPIXj);

  Eigen::MatrixXd calc_wgt_image(const int &NPIXi, const int &NPIXj)
  {
    Eigen::MatrixXd result(Eigen::MatrixXd::Zero(NPIXi,NPIXj));
    for(size_t r=0;r<responses.size();++r)
      {
        for(int j=j0_im[r];j<j1_im[r];++j)
          for(int i=i0_im[r];i<i1_im[r];++i)
            {
              result(i,j)+=responses[r](i+i0_ft[r]-i0_im[r],
                                        j+j0_ft[r]-j0_im[r]);
            }
      }
    LOG4CXX_INFO(logger,"Weight array computed\n");
    return result;
  }
};

#endif
