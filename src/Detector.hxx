#pragma once

#include <string>
#include <valarray>
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <CCfits/CCfits>
#include "Exception.hxx"

namespace hires
{
class Detector
{
public:
  int id;
  double radius_radians, pix_per_radian;
  int nx, ny;
  std::valarray<double> detector_response;

  Detector () {}
  Detector (const boost::filesystem::path &p)
  {
    CCfits::FITS drf (p.string (), CCfits::Read, true);
    CCfits::PHDU &phdu (drf.pHDU ());
    phdu.readAllKeys ();
    nx = phdu.axis (0);
    ny = phdu.axis (1);
    std::map<std::string, CCfits::Keyword *> &m (phdu.keyWord ());
    m["DETECTOR"]->value (id);
    double cdelt1, cdelt2;
    m["CDELT1"]->value (cdelt1);
    m["CDELT2"]->value (cdelt2);

    if (nx != ny)
      throw Exception ("In " + p.string () + " NAXIS1 must equal NAXIS2\n");
    if (nx % 2 == 0)
      throw Exception ("In " + p.string () + " NAXIS1 must be odd\n");
    // FIXME: This seems bogus.  Floating point comparisons can be tricky.
    if (std::abs (cdelt1) != std::abs (cdelt2))
      throw Exception ("In " + p.string ()
                       + " CDELT1 and CDELT2 must be same size\n");
    if (cdelt1 >= 0)
      throw Exception ("In " + p.string () + " CDELT1 must be negative\n");
    if (cdelt2 <= 0)
      throw Exception ("In " + p.string () + " CDELT2 must be positive\n");

    pix_per_radian = 1/(cdelt2 * boost::math::constants::pi<double>() / 180);
    int radius_pix = nx / 2;
    radius_radians = radius_pix / pix_per_radian;

    phdu.read (detector_response);
  }

  double response (const double &du, const double &dv) const
  {
    double result (0);
    int i ((du + radius_radians) * pix_per_radian + 0.5),
        j ((dv + radius_radians) * pix_per_radian + 0.5);
    if (!(i < 0 || i >= nx || j < 0 || j >= ny))
      {
        // FIXME: Column or row major?
        result = detector_response[i + nx * j];
      }
    return result;
  }
};
}
