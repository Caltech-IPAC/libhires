#ifndef HIRES_FOOTPRINT_HXX
#define HIRES_FOOTPRINT_HXX

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

  std::vector<double> flux;
  std::vector<int> j0_im, j1_im, i0_im, i1_im, j0_ft, j1_ft, i0_ft, i1_ft;

  Footprint (const double &radians_per_pix, const int &ni, const int &nj,
             const double &min_sample_flux, const double &angle_tolerance,
             const double &footprints_per_pix,
             const std::map<int, Detector> &detectors,
             std::vector<Sample> &samples);

  double count_good_samples (const double &radians_per_pix, const int &ni,
                             const int &nj,
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
                                   const int &ni, const int &nj) const;

  void compute_correction (const int &nx, const int &ny,
                           const arma::mat &flux_image, const int &iter,
                           const bool &do_cfv,
                           const std::function<double(double)> &boost_func,
                           const int &boost_max_iter, arma::mat &correction,
                           arma::mat &correction_squared) const;

  void compute_minimap (const double &radians_per_pix, const int &nx,
                        const int &ny, const std::vector<Sample> &samples,
                        arma::mat &minimap_image, arma::mat &minimap_hitmap) const;

  void set_fluxes_to_sim_values (const arma::mat &sim_image);

  arma::mat calc_wgt_image (const int &ni, const int &nj) const
  {
    arma::mat result (nj, ni);
    result.fill (1e-8);
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

#endif
