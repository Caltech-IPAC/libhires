#include "../Hires.hxx"
#include "../Footprint.hxx"

namespace hires
{
Hires::Hires (const std::array<size_t,2> &Nxy,
              const std::array<double,2> &Crval, const double &Radians_per_pix,
              const std::vector<std::pair<std::string, std::pair<std::string,
                                                                 std::string> > >
              &Fits_keywords,
              const std::vector<Sample> &Samples):
  nxy(Nxy),
  crval(Crval), radians_per_pix(Radians_per_pix),
  angle_tolerance (2.5),
  fits_keywords(Fits_keywords),
  samples(Samples)
{
  fits_keywords.push_back({std::string ("CREATOR"),
        {std::string ("LIBHIRES"), std::string("")}});
}
}
