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
  Eigen::MatrixXd response;

  std::vector<double> signal;
  std::vector<int> j0_im, j1_im, i0_im, i1_im, j0_ft, j1_ft, i0_ft, i1_ft;
  Eigen::VectorXd result;
  Eigen::VectorXd rhs;

  double noise_level;
  double lambda;

  Footprint (const double &radians_per_pix, const std::array<int,2> &nxy,
             const double &angle_tolerance,
             const double &footprints_per_pix,
             const std::map<int, Detector> &detectors,
             const std::vector<Sample> &samples);
  Footprint () {}

  double count_good_samples (const double &radians_per_pix,
                             const std::array<int,2> &nxy,
                             const std::vector<Sample> &samples);
  Eigen::MatrixXd get_response (const int &detector_id, const double &i_frac,
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

  void compute_correction (const Eigen::MatrixXd &signal_image,
                           Eigen::MatrixXd &correction);

};
}
