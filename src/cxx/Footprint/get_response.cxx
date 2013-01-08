#include <tuple>
#include <vector>
#include "../Footprint.hxx"

/* Return appropriate footprint array, generating it if needed */

std::vector<double>
Footprint::get_response(const int &detector_id, const double &i_frac,
                        const double &j_frac, const double &angle,
                        const double &angle_tolerance,
                        const double &footprints_per_pix,
                        const std::map<int,Detector> &detectors)
{
  double rounded_angle=std::round(angle/angle_tolerance);
  const int angle_id(rounded_angle);
  double recomposed_angle(rounded_angle*angle_tolerance);
  int i_id(i_frac*footprints_per_pix), j_id(j_frac*footprints_per_pix);
  auto key=std::make_tuple(detector_id,angle_id,i_id,j_id);
  auto iter(responses.find(key));
  if(iter==responses.end())
    {
      double delta(1.0/(footprints_per_pix)), zero(delta/2.0 -0.5);
      double i_mod(zero+i_id*delta), j_mod(zero+j_id*delta);
      /* Ignore the result, since we know that the element is not
         already in the map */
      iter=responses.insert
        (std::make_pair(key,generate_response
                        (detector_id,i_mod,j_mod,recomposed_angle,
                         detectors))).first;
    }
  return iter->second;
}
