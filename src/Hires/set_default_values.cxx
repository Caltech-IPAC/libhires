#include <boost/math/constants/constants.hpp>

#include "../Hires.hxx"

void hires::Hires::set_default_values ()
{
  drf_prefix = "";
  ctype1 = "RA -- TAN";
  ctype2 = "DEC -- TAN";
  boost_type = "TIMES_2";
  flux_units = "??";
  outfile_types = { "hires" }, ni = 500;
  nj = 500;
  boost_max_iter = 0;
  footprints_per_pix = 1;
  beam_spike_n = 5;
  radians_per_pix = 60 * boost::math::constants::pi<double>() / (3600 * 180);
  crval1 = std::numeric_limits<double>::max();
  crval2 = std::numeric_limits<double>::max();
  min_sample_flux = std::numeric_limits<double>::min ();
  angle_tolerance = 2.5;
  beam_spike_height = 10;

  std::string kwd ("AUTHOR");
  std::string val ("LIBHIRES");
  std::string comment ("");
  fits_keywords.push_back (std::make_tuple (kwd, val, comment));

  boost_func = [](const double &x)
  { return x + x - 1.0; };
}