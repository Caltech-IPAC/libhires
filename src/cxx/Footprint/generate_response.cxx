/* generate full footprint array for given detector at given angle and
   centered at given x,y offset from pixel center.
   -.5 <= i_offset,j_offset < 0.5 */

#include <vector>
#include "../Footprint.hxx"

Eigen::MatrixXd
Footprint::generate_response(const int &detector_id, const double &i_offset,
                             const double &j_offset,
                             const double &recomposed_angle,
                             const std::map<int,Detector> &detectors)
{
  const Detector &detector(detectors.find(detector_id)->second);
  const int radius_pix(detector.radius_radians/detector.radians_per_pix);
  const double du_offset(i_offset*detector.radians_per_pix),
    dv_offset(j_offset*detector.radians_per_pix);

  int n_ij=2*radius_pix + 1;
  const double radius=radius_pix*detector.radians_per_pix;

  Eigen::MatrixXd response(n_ij,n_ij);
  const double cos_angle(std::cos(recomposed_angle)),
    sin_angle(std::sin(recomposed_angle));
  double sum(0), dx(-radius);
  for(int i=0;i<n_ij;++i)
    {
      double dy(-radius);
      for(int j=0;j<n_ij;++j)
        {
          response(i,j)=
            detector.response(dx*cos_angle - dy*sin_angle - du_offset,
                              dy*cos_angle + dx*sin_angle - dv_offset);
          sum+=response(i,j);
          dy+=detector.radians_per_pix;
        }
      dx+=detector.radians_per_pix;
    }
  response/=sum;
  return response;
}
