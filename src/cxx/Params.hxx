#ifndef HIRES_PARAMS_HXX
#define HIRES_PARAMS_HXX

#include <string>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <functional>

class Params
{
public:
  std::string infile_prefix, outfile_prefix, drf_prefix, ctype1, ctype2,
    boost_type, starting_image, beam_starting_image, flux_units, log_filename;
  std::vector<std::string> outfile_types;
  int NPIXi, NPIXj, iter_max, boost_max_iter, footprints_per_pix, beam_spike_n;
  double radians_per_pix, crval1, crval2, min_sample_flux, angle_tolerance,
    beam_spike_height;
  std::vector<int> iter_list;
  std::vector<std::tuple<std::string,std::string,std::string>> fits_keywords;
  std::function<double (double)> boost_func;

  Params(int argc, char* argv[]);  

  enum class Data_Type {planck, spire} data_type;

  static std::string type_string(const Data_Type &dt)
  {
    if(dt==Data_Type::planck)
      return "planck";
    else if(dt==Data_Type::spire)
      return "spire";
    abort();
    return "";
  }
};

std::ostream& operator<<(std::ostream& out, const Params &p);

#endif
