#ifndef HIRES_PARAMS_HXX
#define HIRES_PARAMS_HXX

#include <string>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <functional>
#include <armadillo>
#include <map>
#include "Hires_Parameters.hxx"
#include "Sample.hxx"
#include "Exception.hxx"

namespace hires
{

class Hires : public Hires_Parameters
{
public:
  arma::mat hitmap, minimap;
  arma::mat wgt_image;
  std::map<int, arma::mat> flux_images;
  std::map<int, arma::mat> cfv_images;
  std::map<int, arma::mat> beam_images;

  // Constructors
  // takes Hires_Parameters as single argument
  Hires (const Hires_Parameters &hp)
      : Hires_Parameters (
            hp.drf_prefix, hp.ctype1, hp.ctype2, hp.boost_type,
            hp.starting_image, hp.beam_starting_image, hp.flux_units,
            hp.log_filename, hp.outfile_types, hp.ni, hp.nj, hp.boost_max_iter,
            hp.footprints_per_pix, hp.beam_spike_n, hp.radians_per_pix,
            hp.crval1, hp.crval2, hp.min_sample_flux, hp.angle_tolerance,
            hp.beam_spike_height, hp.fits_keywords, hp.boost_func) {};

  void iterate (int &iter, std::vector<hires::Sample> &samples);

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
#endif
