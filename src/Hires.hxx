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

namespace hires
{
class Hires
{
public:
  std::array<size_t,2> nxy;
  std::array<double,2> crval;
  double radians_per_pix, angle_tolerance;

  std::vector<std::pair<std::string, std::pair<std::string, std::string> > >
  fits_keywords;
  boost::filesystem::path drf_file;
  
  const std::vector<Sample> &samples;

  arma::mat minimap, hires, elastic_net;

  Hires (const std::array<size_t,2> &Nxy,
         const std::array<double,2> &Crval, const double &Radians_per_pix,
         const std::vector<std::pair<std::string, std::pair<std::string,
                                                            std::string> > >
         &Fits_keywords,
         const std::vector<Sample> &Samples);

  void dump_params ();
  void compute_hires (const boost::filesystem::path &Drf_file);
  void compute_minimap ();

  void write_minimap (const std::string &outfile_prefix)
  {
    write_file(outfile_prefix,"minimap","Minimap Image",false,minimap);
  }
  void write_hires (const std::string &outfile_prefix)
  {
    write_file(outfile_prefix,"hires","HIRES image",true,hires);
  }
  void write_file (const std::string &output_prefix,
                   const std::string &filename,
                   const std::string &filetype,
                   const bool add_drf_filename,
                   const arma::mat &image);
                   
  void write_fits (const arma::mat &image,
                   const std::vector<std::pair<std::string,
                                               std::pair<std::string,
                                                         std::string> > >
                   &file_specific_keywords, const std::string &outfile_name);
};

std::ostream &operator<<(std::ostream &out, const Hires &p);
}
