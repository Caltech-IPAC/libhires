#pragma once

#include <string>
#include <iostream>
#include <functional>
#include <map>

#include <armadillo>
#include <boost/algorithm/string.hpp>
#include <boost/math/constants/constants.hpp>

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

  size_t iteration;

  arma::mat hitmap, minimap;
  arma::mat wgt_image;
  std::map<int, arma::mat> flux_images;
  std::map<int, arma::mat> cfv_images;
  std::map<int, arma::mat> beam_images;

  Hires (): ctype1 ("RA -- TAN"), ctype2 ("DEC -- TAN"),
            boost_type ("TIMES_2"), flux_units ("??"),
            outfile_types ({ "hires" }), ni (500),
            nj (500), boost_max_iter (0),
            footprints_per_pix (1), beam_spike_n (5),
            radians_per_pix (60 * boost::math::constants::pi<double>() / (3600 * 180)),
            crval1 (std::numeric_limits<double>::max()),
            crval2 (std::numeric_limits<double>::max()),
            min_sample_flux (std::numeric_limits<double>::min ()),
            angle_tolerance (2.5),
            beam_spike_height (10),
            fits_keywords ({std::make_tuple
                  (std::string ("AUTHOR"),
                   std::string ("LIBHIRES"),
                   std::string(""))}),
            boost_func ([](const double &x)
                        { return x + x - 1.0; }),
            iteration(0) {}

  Hires (const std::vector<std::string> param_str) : Hires ()
  {
    parse_command_line (param_str);
  }


  void parse_command_line (const std::vector<std::string> param_str);
  void dump_params ();

  void init (std::vector<hires::Sample> &samples);
  void iterate (std::vector<hires::Sample> &samples);

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

  void write_output (const Image_Type image_type,
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
