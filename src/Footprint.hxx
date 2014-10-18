#pragma once

#include <array>
#include <vector>
#include <valarray>
#include <map>
#include <armadillo>
#include "Detector.hxx"
#include "Sample.hxx"

namespace hires
{
class Footprint
{
public:
  std::vector<std::valarray<bool> > good;
  std::map<std::tuple<int, int, int, int>, arma::mat> responses_complete;
  std::vector<const arma::mat *> responses;

  std::vector<double> signal;
  std::vector<int> j0_im, j1_im, i0_im, i1_im, j0_ft, j1_ft, i0_ft, i1_ft;

  Footprint (const double &radians_per_pix, const std::array<int,2> &nxy,
             const double &min_sample_signal, const double &angle_tolerance,
             const double &footprints_per_pix,
             const std::map<int, Detector> &detectors,
             std::vector<Sample> &samples);

  double count_good_samples (const double &radians_per_pix,
                             const std::array<int,2> &nxy,
                             const std::vector<Sample> &samples);
  const arma::mat *get_response (const int &detector_id, const double &i_frac,
                                 const double &j_frac, const double &angle,
                                 const double &angle_tolerance,
                                 const double &footprints_per_pix,
                                 const double &radians_per_pix,
                                 const std::map<int, Detector> &detectors);

  arma::mat generate_response (const int &detector_id, const double &i_offset,
                               const double &j_offset,
                               const double &recomposed_angle,
                               const double &radians_per_pix,
                               const std::map<int, Detector> &detectors) const;

  std::vector<int> compute_bounds (const arma::mat &response,
                                   const int &i_center, const int &j_center,
                                   const std::array<int,2> &nxy) const;

  void compute_correction (const std::array<int,2> &nxy,
                           const arma::mat &signal_image, const int &iter,
                           const bool &do_cfv, const bool &boosting,
                           const std::function<double(double)> &boost_function,
                           arma::mat &correction,
                           arma::mat &correction_squared) const;

  void compute_minimap (const double &radians_per_pix,
                        const std::array<int,2> &nxy,
                        const std::vector<Sample> &samples,
                        arma::mat &minimap_image, arma::mat &minimap_hitmap) const;

  void set_signals_to_sim_values (const arma::mat &sim_image);

  arma::mat calc_wgt_image (const std::array<int,2> &nxy) const
  {
    arma::mat result (nxy[1], nxy[0]);
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
