#pragma once

#include <array>
#include <string>
#include <iostream>
#include <functional>
#include <map>
#include <set>

#include <armadillo>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

#include "Sample.hxx"
#include "Exception.hxx"
#include "Detector.hxx"
#include "Footprint.hxx"

namespace hires
{
std::map<int, Detector> read_DRF (const boost::filesystem::path &DRF_file);

class Hires
{
public:
  enum class Image_Type
  {
    hires_image,
    hires_covariance,
    hires_correction,
    hires_beam,
    minimap_image,
    minimap_hitmap
  };

  std::string starting_image, beam_starting_image;
  boost::filesystem::path drf_file;
  std::array<int,2> nxy;
  int footprints_per_pix, beam_spike_n;
  std::array<double,2> crval;
  double radians_per_pix, angle_tolerance, beam_spike_height;
  std::set<Image_Type> output_types;

  std::vector<std::pair<std::string, std::pair<std::string, std::string> > >
  fits_keywords;

  std::function<double(double)> boost_function;
  
  // FIXME: Why isn't this const?
  const std::vector<Sample> &samples;

  std::map<int, Detector> detectors;
  Footprint footprints;
  size_t iteration;
  arma::mat hitmap, minimap, weight_image, signal_image;

  Hires (const std::array<int,2> &Nxy,
         const std::array<double,2> &Crval, const double &Radians_per_pix,
         const std::set<Image_Type> &Output_types,
         const boost::filesystem::path &Drf_file,
         const std::string &boost_function_string,
         const std::vector<std::pair<std::string, std::pair<std::string,
                                                            std::string> > >
         &Fits_keywords,
         const std::vector<Sample> &Samples);

  void dump_params ();
  void iterate (const bool &boosting);
  void compute_minimap ();
  bool running_hires() const
  {
    return output_types.find(Image_Type::hires_image)!=output_types.end()
      || output_types.find(Image_Type::hires_covariance)!=output_types.end()
      || output_types.find(Image_Type::hires_correction)!=output_types.end()
      || output_types.find(Image_Type::hires_beam)!=output_types.end();
  }

  void write_output (const std::string &outfile_prefix);
  void write_file (const std::string &output_prefix, const Image_Type &type);

  arma::mat start_image (const std::string &filename, int &iter_start);

  void write_fits (const arma::mat &image,
                   const std::vector<std::pair<std::string,
                                               std::pair<std::string,
                                                         std::string> > >
                   &file_specific_keywords, const std::string &outfile_name);
};

std::ostream &operator<<(std::ostream &out, const Hires &p);
}
