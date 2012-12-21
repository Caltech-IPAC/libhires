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
  double deg_per_pix, crval1, crval2, min_sample_flux, angle_tolerance,
    beam_spike_height;
  std::vector<int> iter_list;
  std::vector<std::tuple<std::string,std::string,std::string>> fits_keywords;
  std::function<double (double)> boost_func;

  Params(int argc, char* argv[]);  

  enum class Data_Type {planck, spire} data_type;
};

std::ostream& operator<<(std::ostream& out, const Params &p);

#endif
