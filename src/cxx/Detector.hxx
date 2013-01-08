#ifndef HIRES_DETECTOR_HXX
#define HIRES_DETECTOR_HXX

#include <string>
#include <valarray>
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <CCfits>
#include "logger.hxx"

class Detector
{
public:
  int id;
  double radius_radians, radians_per_pix;
  int nx, ny;
  std::valarray<double> detector_response;

  Detector() {}
  Detector(const boost::filesystem::path &p)
  {
    CCfits::FITS drf(p.string(),CCfits::Read,true);
    CCfits::PHDU &phdu(drf.pHDU());
    phdu.readAllKeys();
    nx=phdu.axis(0);
    ny=phdu.axis(1);
    std::map<std::string,CCfits::Keyword *> &m(phdu.keyWord());
    m["DETECTOR"]->value(id);
    double cdelt1, cdelt2;
    m["CDELT1"]->value(cdelt1);
    m["CDELT2"]->value(cdelt2);

    if(nx!=ny)
      LOG4CXX_ERROR(logger,"In " << p.string() << " NAXIS1 must equal NAXIS2\n");
    if(nx%2==0)
      LOG4CXX_ERROR(logger,"In " << p.string() << " NAXIS1 must be odd\n");
    /* This seems bogus.  Floating point comparisons can be tricky. */
    if(abs(cdelt1) != abs(cdelt2))
      LOG4CXX_ERROR(logger, "In " << p.string()
                    << " CDELT1 and CDELT2 must be same size\n");
    if(cdelt1 >= 0)
      LOG4CXX_WARN(logger, "In " << p.string() << " CDELT1 must be negative\n");
    if(cdelt2 <= 0)
      LOG4CXX_WARN(logger, "In " << p.string() << " CDELT2 must be positive\n");

    double radians_per_pix=cdelt2*boost::math::constants::pi<double>()/180;
    int radius_pix = nx / 2;
    double radius_radians = radius_pix * radians_per_pix;

    phdu.read(detector_response);
    LOG4CXX_INFO(logger, "detector: " << id 
                 << "; file=" << p.filename()
                 << "; radius= " << radius_pix
                 << "; pixels = "
                 << radius_radians*60*180/boost::math::constants::pi<double>()
                 << " arcmin\n");
  }

  double response(const double &du, const double &dv) const
  {
    double result(0);
    int i((du+radius_radians)/radians_per_pix + 0.5),
      j((dv+radius_radians)/radians_per_pix + 0.5);
    if(!(i<0 || i>=nx || j<0 || j>ny))
      {
        /* Column or row major? */
        result=detector_response[i+nx*j];
      }
    return result;
  }
};

namespace std
{
  template<>
  inline void swap(Detector &a, Detector &b)
  {
    swap(a.id,b.id);
    swap(a.radius_radians,b.radius_radians);
    swap(a.radians_per_pix,b.radians_per_pix);
    swap(a.nx,b.nx);
    swap(a.ny,b.ny);
    swap(a.detector_response,b.detector_response);
  }
}

#endif