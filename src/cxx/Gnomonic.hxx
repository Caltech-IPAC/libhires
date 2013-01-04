#ifndef HIRES_GNOMONIC_HXX
#define HIRES_GNOMONIC_HXX

/* Compute gnomonic (tangent plane) projection of lon,lat to x,y
   lat,lon,x,y all in degrees
   Equations from: http://mathworld.wolfram.com/GnomonicProjection.html
   Example:
   gn = Gnomonic(lon_center_degrees, lat_center_degrees)
   x, y = lonlat2xy(lon_degrees, lat_degrees)
*/

#include <cmath>
#include <boost/math/constants/constants.hpp>

class Gnomonic
{
public:
  double lon0, sin_lat0, cos_lat0;

  Gnomonic(const double &lon_degrees, const double &lat_degrees):
    lon0(radians(lon_degrees)),
    sin_lat0(std::sin(radians(lat_degrees))),
    cos_lat0(std::cos(radians(lat_degrees))) {}
  Gnomonic()=delete;

  double radians(const double &degree) const
  {
    return degree*boost::math::constants::pi<double>()/180;
  }

  template<class T>
  std::pair<T,T> lonlat2xy(const T &lon, const T &lat) const
  {
    T cos_dlon(std::cos(lon - lon0)), sin_dlon(std::sin(lon - lon0));
    T sin_lat(std::sin(lat)), cos_lat(std::cos(lat));
    T cos_c(sin_lat0*sin_lat + cos_lat0*cos_lat*cos_dlon);
    T x = (cos_lat * sin_dlon) / cos_c;
    T y = (cos_lat0 * sin_lat - sin_lat0 * cos_lat * cos_dlon) / cos_c;
    return std::make_pair(x, y);
  }
};

#endif

