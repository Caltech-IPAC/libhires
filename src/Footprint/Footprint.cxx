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

  std::array<double,2> offset{{nxy[0] / 2.0, nxy[1] / 2.0}};

  for (size_t s = 0; s < samples.size (); ++s)
    {
      const size_t n (good[s].size ());
      for (size_t j = 0; j < n; ++j)
        {
          if (!good[s][j])
            continue;
          double xi ((samples[s].x[j] / radians_per_pix) + offset[0]),
              yi ((samples[s].y[j] / radians_per_pix) + offset[1]);

          int i_int (xi), j_int (yi);
          double i_frac (xi - i_int), j_frac (yi - j_int);

          responses.push_back (
              get_response (samples[s].id, i_frac, j_frac, samples[s].angle[j],
                            angle_tolerance, footprints_per_pix,
                            radians_per_pix, detectors));
          std::vector<int> bounds (
              compute_bounds (*(*responses.rbegin ()), i_int, j_int, nxy));

          signal.push_back (samples[s].signal[j]);
          j0_im.push_back (bounds[0]);
          j1_im.push_back (bounds[1]);
          i0_im.push_back (bounds[2]);
          i1_im.push_back (bounds[3]);
          j0_ft.push_back (bounds[4]);
          j1_ft.push_back (bounds[5]);
          i0_ft.push_back (bounds[6]);
          i1_ft.push_back (bounds[7]);
        }
    }

  const size_t num_pixels=nxy[0]*nxy[1];
  Eigen::MatrixXd response(num_pixels,num_pixels);
  response.setZero();
  Eigen::VectorXd rhs(num_pixels);
  rhs.setZero();

  std::cout << "num_pixels: " << nxy[0] << " "
            << nxy[1] << " "
            << num_pixels << " "
            << signal.size() << " "
            << "\n";
  for (size_t n = 0; n < signal.size(); ++n)
    {

          if(n%10000==0)
            std::cout << "signal: " << n << "\n";

      for (int column_i = i0_im[n]; column_i < i1_im[n] ; ++column_i)
           
        for (int column_j = j0_im[n]; column_j < j1_im[n]; ++column_j)
          {
            double column_response=
              (*responses[n])(column_j - j0_im[n] + j0_ft[n],
                              column_i - i0_im[n] + i0_ft[n]);
            const size_t column=column_i + nxy[0]*column_j;
            rhs(column)=signal[n]*column_response;
            for (int row_j = j0_im[n]; row_j < j1_im[n]; ++row_j)
              for (int row_i = i0_im[n]; row_i < i1_im[n]; ++row_i)
                {
                  const size_t row=row_i + nxy[0]*row_j;
                  response(row,column)+=column_response
                    * (*responses[n])(row_j - j0_im[n] + j0_ft[n],
                                      row_i - i0_im[n] + i0_ft[n]);
                }
          }
    }
  std::cout << "solving Svd\n";

  std::cout << "response:\n"
            << response.determinant() << "\n"
            << (response.transpose() * response).determinant() << "\n";
            // << rhs << "\n";
  // result=response.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
  // std::cout << "solving QR\n";
  // result=response.householderQr().solve(rhs);
  // std::cout << "solving colQR\n";
  // result=response.colPivHouseholderQr().solve(rhs);
  // std::cout << "solving fullQR\n";
  // result=response.fullPivHouseholderQr().solve(rhs);
  std::cout << "solving Normal\n";
  result=(response.transpose() * response).ldlt().solve(response.transpose()*rhs);
  std::cout << "solved\n";
  // std::cout << "result:\n"
  //           << result << "\n";

}
}
