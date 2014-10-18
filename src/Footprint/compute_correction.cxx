#include <armadillo>
#include "../Footprint.hxx"

namespace hires
{
void Footprint::compute_correction (
    const std::array<int,2> &nxy, const arma::mat &signal_image, const int &iter,
    const bool &do_cfv, const bool &boosting, 
    const std::function<double(double)> &boost_function,
    arma::mat &correction,
    arma::mat &correction_squared) const
{
  correction.zeros (nxy[1], nxy[0]);
  if (do_cfv)
    correction_squared.zeros (nxy[1], nxy[0]);

  for (size_t n = 0; n < signal.size (); ++n)
    {
      arma::mat integration (j1_im[n] - j0_im[n], i1_im[n] - i0_im[n]);
      for (int i = 0; i < i1_im[n] - i0_im[n]; ++i)
        for (int j = 0; j < j1_im[n] - j0_im[n]; ++j)
          {
            integration (j, i) = signal_image (j + j0_im[n], i + i0_im[n])
                                 * (*responses[n])(j + j0_ft[n], i + i0_ft[n]);
          }

      double signal_prime (accu (integration));

      double scale (signal[n] / signal_prime);
      if (iter != 1 && boosting && boost_function)
        {
          scale = boost_function (scale);
        }
      for (int i = 0; i < i1_im[n] - i0_im[n]; ++i)
        for (int j = 0; j < j1_im[n] - j0_im[n]; ++j)
          {
            correction (j + j0_im[n], i + i0_im[n])
                += (*responses[n])(j + j0_ft[n], i + i0_ft[n]) * scale;
          }

      if (do_cfv)
        {
          for (int i = 0; i < i1_im[n] - i0_im[n]; ++i)
            for (int j = 0; j < j1_im[n] - j0_im[n]; ++j)
              {
                correction_squared (j + j0_im[n], i + i0_im[n])
                    += (*responses[n])(j + j0_ft[n], i + i0_ft[n])
                       * scale *scale;
              }
        }
    }
}
}
