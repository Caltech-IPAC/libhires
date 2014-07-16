#pragma once

#include <string>
#include <iostream>
#include <functional>
#include <map>

#include <armadillo>
#include <boost/algorithm/string.hpp>

#include "Sample.hxx"
#include "Exception.hxx"

namespace hires
{

class Hires
{
public:
  std::string drf_prefix, ctype1, ctype2, boost_type, starting_image,
      beam_starting_image, flux_units;
  std::vector<std::string> outfile_types;
  int ni, nj, boost_max_iter, footprints_per_pix, beam_spike_n;
  double radians_per_pix, crval1, crval2, min_sample_flux, angle_tolerance,
      beam_spike_height;

  std::vector<std::tuple<std::string, std::string, std::string> >
  fits_keywords;
  std::function<double(double)> boost_func;

  arma::mat hitmap, minimap;
  arma::mat wgt_image;
  std::map<int, arma::mat> flux_images;
  std::map<int, arma::mat> cfv_images;
  std::map<int, arma::mat> beam_images;

  Hires (std::string drf_prefix, std::string ctype1, std::string ctype2,
         std::string boost_type, std::string starting_image,
         std::string beam_starting_image, std::string flux_units,
         std::vector<std::string> outfile_types,
         int ni, int nj, int boost_max_iter, int footprints_per_pix,
         int beam_spike_n, double radians_per_pix, double crval1,
         double crval2, double min_sample_flux, double angle_tolerance,
         double beam_spike_height,
         std::vector<std::tuple<std::string, std::string, std::string> >
             fits_keywords, std::function<double(double)> boost_func)
      : drf_prefix (drf_prefix), ctype1 (ctype1), ctype2 (ctype2),
        boost_type (boost_type), starting_image (starting_image),
        beam_starting_image (beam_starting_image), flux_units (flux_units),
        outfile_types (outfile_types), ni (ni),
        nj (nj), boost_max_iter (boost_max_iter),
        footprints_per_pix (footprints_per_pix), beam_spike_n (beam_spike_n),
        radians_per_pix (radians_per_pix), crval1 (crval1), crval2 (crval2),
        min_sample_flux (min_sample_flux), angle_tolerance (angle_tolerance),
        beam_spike_height (beam_spike_height), fits_keywords (fits_keywords),
        boost_func (boost_func)
  {
  }

  Hires (const std::vector<std::string> param_str) : Hires ()
  {
    parse_command_line (param_str);
  }

  Hires () { set_default_values (); }

  void set_default_values ();
  void parse_command_line (const std::vector<std::string> param_str);
  void dump_params ();

  void iterate (int &iter, std::vector<hires::Sample> &samples);

  enum class Image_Type
  {
    all,
    hires,
    cov,
    cfv,
    beam,
    minimap,
    hitmap
  };

  void write_output (int &iter, const Image_Type image_type,
                     const std::string &outfile_prefix);

  void write_file (arma::mat image, std::string filename, const char *desc,
                   int iter, int isflux);

  arma::mat spike_image ();
  arma::mat start_image (const std::string &filename, int &iter_start);

  void write_fits (
      const arma::mat &image,
      const std::vector<std::tuple<std::string, std::string, std::string> >
          file_specific_keywords, const std::string &outfile_name);
};

std::ostream &operator<<(std::ostream &out, const Hires &p);
}
