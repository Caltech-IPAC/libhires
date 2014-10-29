#include <cmath>

#include "../Footprint.hxx"

namespace hires
{
Footprint::Footprint (const double &radians_per_pix,
                      const std::array<int,2> &nxy,
                      const double &angle_tolerance,
                      const double &footprints_per_pix,
                      const std::map<int, Detector> &detectors,
                      const std::vector<Sample> &samples)
{
  size_t num_good = count_good_samples (radians_per_pix, nxy, samples);
  j0_im.reserve (num_good);
  j1_im.reserve (num_good);
  i0_im.reserve (num_good);
  i1_im.reserve (num_good);
  j0_ft.reserve (num_good);
  j1_ft.reserve (num_good);
  i0_ft.reserve (num_good);
  i1_ft.reserve (num_good);

  const size_t num_pixels=nxy[0]*nxy[1];
  response.resize(num_pixels,num_pixels);
  response.setZero();
  rhs.resize(num_pixels);
  rhs.setZero();

  std::array<double,2> offset{{nxy[0] / 2.0, nxy[1] / 2.0}};

  // FIXME: Get rid of vector of samples and just have a Sample?
  for (size_t s = 0; s < samples.size (); ++s)
    {
      // FIXME: Parallelize this
      for (size_t j = 0; j < good[s].size (); ++j)
        {
          if (!good[s][j])
            continue;
          double xi ((samples[s].x[j] / radians_per_pix) + offset[0]),
              yi ((samples[s].y[j] / radians_per_pix) + offset[1]);

          int i_int (xi), j_int (yi);
          double i_frac (xi - i_int), j_frac (yi - j_int);

          Eigen::MatrixXd detector_response=
            get_response (samples[s].id, i_frac, j_frac, samples[s].angle[j],
                          angle_tolerance, footprints_per_pix,
                          radians_per_pix, detectors);

          std::vector<int> bounds (compute_bounds (detector_response, i_int,
                                                   j_int, nxy));

          signal.push_back (samples[s].signal[j]);
          j0_im.push_back (bounds[0]);
          j1_im.push_back (bounds[1]);
          i0_im.push_back (bounds[2]);
          i1_im.push_back (bounds[3]);
          j0_ft.push_back (bounds[4]);
          j1_ft.push_back (bounds[5]);
          i0_ft.push_back (bounds[6]);
          i1_ft.push_back (bounds[7]);

          /// Renormalize the response function by what intersects
          /// with the image

          double sum=0;
          for(int iy=*(j0_ft.rbegin()); iy<*(j1_ft.rbegin()); ++iy)
            for(int ix=*(i0_ft.rbegin()); ix<*(i1_ft.rbegin()); ++ix)
              sum+=detector_response(ix,iy);
          detector_response/=sum;

          for (int column_i = *i0_im.rbegin(); column_i < *i1_im.rbegin();
               ++column_i)
            for (int column_j = *j0_im.rbegin(); column_j < *j1_im.rbegin();
                 ++column_j)
              {
                double column_response=
                  detector_response(column_j - *j0_im.rbegin() + *j0_ft.rbegin(),
                               column_i - *i0_im.rbegin() + *i0_ft.rbegin());
                const size_t column=column_i + nxy[0]*column_j;
                rhs(column)+=*signal.rbegin()*column_response;

                for (int row_j = *j0_im.rbegin(); row_j < *j1_im.rbegin();
                     ++row_j)
                  for (int row_i = *i0_im.rbegin(); row_i < *i1_im.rbegin();
                       ++row_i)
                    {
                      const size_t row=row_i + nxy[0]*row_j;
                      response(row,column)+=column_response
                        * detector_response(row_j - *j0_im.rbegin() + *j0_ft.rbegin(),
                                            row_i - *i0_im.rbegin() + *i0_ft.rbegin());
                    }
              }
        }
    }

  result.resize(nxy[0]*nxy[1]);
  result.setZero();
  
  // FIXME: Making this number smaller seems to break things
  noise_level=1e-2;

  // FIXME: I got this from dimensional analysis and added the 0.5
  // because I felt like it.
  lambda=0.5*noise_level*noise_level*signal.size()/response.rows();
}
}
