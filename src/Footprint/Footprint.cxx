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

          Eigen::MatrixXd response=get_response (samples[s].id, i_frac, j_frac, samples[s].angle[j],
                                                 angle_tolerance, footprints_per_pix,
                                                 radians_per_pix, detectors);

          std::vector<int> bounds (compute_bounds (response, i_int, j_int,
                                                   nxy));

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
              sum+=response(ix,iy);
          response/=sum;

          // FIXME: Use std::move?
          responses.push_back (response);
        }
    }

  const size_t num_pixels=nxy[0]*nxy[1];
  Eigen::MatrixXd response(num_pixels,num_pixels);
  response.setZero();
  Eigen::VectorXd rhs(num_pixels);
  rhs.setZero();

  // FIXME: Parallelize this
  for (size_t n = 0; n < signal.size(); ++n)
    {

          if(n%10000==0)
            std::cout << "signal: " << n << "\n";

      for (int column_i = i0_im[n]; column_i < i1_im[n] ; ++column_i)
           
        for (int column_j = j0_im[n]; column_j < j1_im[n]; ++column_j)
          {
            double column_response=
              responses[n](column_j - j0_im[n] + j0_ft[n],
                           column_i - i0_im[n] + i0_ft[n]);
            const size_t column=column_i + nxy[0]*column_j;
            rhs(column)+=signal[n]*column_response;

            for (int row_j = j0_im[n]; row_j < j1_im[n]; ++row_j)
              for (int row_i = i0_im[n]; row_i < i1_im[n]; ++row_i)
                {
                  const size_t row=row_i + nxy[0]*row_j;
                  response(row,column)+=column_response
                    * responses[n](row_j - j0_im[n] + j0_ft[n],
                                   row_i - i0_im[n] + i0_ft[n]);
                }
          }
    }

  result.resize(nxy[0]*nxy[1]);
  result.setZero();
  
  // FIXME: Making this number smaller seems to break things
  const double noise_level=1e-2;

  // FIXME: This needs to be set by the expected noise level
  const double lambda=0.5*noise_level*noise_level*signal.size()/response.rows();
  Eigen::MatrixXd &AA=response;

  Eigen::MatrixXd ddB(result.size(),result.size());
  ddB.setZero();
  Eigen::VectorXd dB(result.size());

  const double epsilon=1e-4;
  double result_norm=0, du_norm=0;
  for(int iteration=0; iteration < 40 && result_norm*epsilon <= du_norm ;
      ++iteration)
    {
      Eigen::VectorXd dA(2*AA*result - rhs);

      Eigen::VectorXd entropy(result.size());

      for(int i=0; i<result.size(); ++i)
        {
          entropy(i)=log(std::cosh(result(i)/noise_level));
          dB(i)=std::tanh(result(i)/noise_level)/noise_level;
          double temp=std::cosh(result(i)/noise_level)*noise_level;
          ddB(i,i)=1/(temp*temp);
        }

      Eigen::VectorXd du=-(AA + lambda*ddB).fullPivHouseholderQr()
        .solve(dA + lambda*dB);

      result+=du;

      result_norm=result.norm();
      du_norm=du.norm();

      std::cout << result_norm << " "
                << du_norm << "\n";
    }
}
}
