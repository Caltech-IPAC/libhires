#pragma once

#include <array>
#include <vector>
#include <valarray>
#include <map>
#include <eigen3/Eigen/Eigen>
#include "Detector.hxx"
#include "Sample.hxx"

namespace hires
{
class Footprint
{
public:
  std::vector<std::valarray<bool> > good;
  std::map<std::tuple<int, int, int, int>, Eigen::MatrixXd> responses_complete;
  std::vector<const Eigen::MatrixXd *> responses;

  std::vector<double> signal;
  std::vector<int> j0_im, j1_im, i0_im, i1_im, j0_ft, j1_ft, i0_ft, i1_ft;

  Footprint (const double &radians_per_pix, const std::array<int,2> &nxy,
             const double &angle_tolerance,
             const double &footprints_per_pix,
             const std::map<int, Detector> &detectors,
             const std::vector<Sample> &samples);
  Footprint () {}

  double count_good_samples (const double &radians_per_pix,
                             const std::array<int,2> &nxy,
                             const std::vector<Sample> &samples);
  const Eigen::MatrixXd *get_response (const int &detector_id, const double &i_frac,
                                 const double &j_frac, const double &angle,
                                 const double &angle_tolerance,
                                 const double &footprints_per_pix,
                                 const double &radians_per_pix,
                                 const std::map<int, Detector> &detectors);

  Eigen::MatrixXd generate_response (const int &detector_id, const double &i_offset,
                               const double &j_offset,
                               const double &recomposed_angle,
                               const double &radians_per_pix,
                               const std::map<int, Detector> &detectors) const;

  std::vector<int> compute_bounds (const Eigen::MatrixXd &response,
                                   const int &i_center, const int &j_center,
                                   const std::array<int,2> &nxy) const;

  void compute_correction (const std::array<int,2> &nxy,
                           const Eigen::MatrixXd &signal_image, const int &iter,
                           const bool &boosting,
                           const std::function<double(double)> &boost_function,
                           Eigen::MatrixXd &correction) const;

  Eigen::MatrixXd calc_wgt_image (const std::array<int,2> &nxy) const
  {
    Eigen::MatrixXd result (nxy[1], nxy[0]);
    const double minimum_weight=1e-8;
    result.fill (minimum_weight);
    for (size_t r = 0; r < responses.size (); ++r)
      {
        for (int i = i0_im[r]; i < i1_im[r]; ++i)
          for (int j = j0_im[r]; j < j1_im[r]; ++j)
            result (j, i) += (*responses[r])(j + j0_ft[r] - j0_im[r],
                                             i + i0_ft[r] - i0_im[r]);
      }
    return result;
  }
};
}
