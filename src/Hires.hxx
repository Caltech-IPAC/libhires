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
#include "Binned_Data.hxx"

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
  
  const std::vector<Sample> &samples;

  arma::mat minimap, hires, elastic_net, tikhonov;
  
  Hires (const std::array<size_t,2> &Nxy,
         const std::array<double,2> &Crval, const double &Radians_per_pix,
         const std::vector<std::pair<std::string, std::pair<std::string,
                                                            std::string> > >
         &Fits_keywords,
         const std::vector<Sample> &Samples);

  void compute_minimap ();
  void compute_mcm (const double &sigma_drf, const size_t &num_iterations);
  void compute_elastic_net (const double &sigma_drf);
  void compute_tikhonov (const double &sigma_drf);
  arma::mat compute_response_function (const double &sigma_drf,
                                       const Binned_Data &binned_data);
  
  void write_minimap (const boost::filesystem::path &outfile_prefix)
  {
    write_file(outfile_prefix,"minimap","Minimap Image",minimap);
  }
  void write_mcm (const boost::filesystem::path &outfile_prefix)
  {
    write_file(outfile_prefix,"mcm","MCM Image",hires);
  }
  void write_elastic_net (const boost::filesystem::path &outfile_prefix)
  {
    write_file(outfile_prefix,"elastic_net","Elastic Net Image",elastic_net);
  }
  void write_tikhonov (const boost::filesystem::path &outfile_prefix)
  {
    write_file(outfile_prefix,"tikhonov","Tikhonov Image",tikhonov);
  }
  void write_file (const boost::filesystem::path &output_prefix,
                   const std::string &filename,
                   const std::string &filetype,
                   const arma::mat &image);
                   
  void write_fits (const arma::mat &image,
                   const std::vector<std::pair<std::string,
                                               std::pair<std::string,
                                                         std::string> > >
                   &file_specific_keywords,
                   const boost::filesystem::path &outfile_name);
};

std::ostream &operator<<(std::ostream &out, const Hires &p);
}
