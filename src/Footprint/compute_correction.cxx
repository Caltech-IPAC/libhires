#include <thread>

#include <armadillo>
#include "../Footprint.hxx"

namespace
{
  void thread_callback(const size_t &i, const size_t &num_threads,
                       const arma::mat &signal_image,
                       const std::vector<const arma::mat *> &responses,
                       const std::vector<double> &signal,
                       const bool &boosting, 
                       const std::function<double(double)> &boost_function,
                       const std::array<int,2> &nxy,
                       const int &iter,
                       const std::vector<int> &j0_im,
                       const std::vector<int> &j1_im,
                       const std::vector<int> &i0_im,
                       const std::vector<int> &i1_im,
                       const std::vector<int> &j0_ft,
                       const std::vector<int> &i0_ft,
                       arma::mat &local_correction)

  {
    local_correction.zeros (nxy[1], nxy[0]);

    const size_t size=signal.size ();
    const size_t begin=(size*i/num_threads), end=(size*(i+1)/num_threads); 

    for (size_t n = begin; n < end; ++n)
      {
        double sum=0;
        for (int i = 0; i < i1_im[n] - i0_im[n]; ++i)
          for (int j = 0; j < j1_im[n] - j0_im[n]; ++j)
            {
              sum += signal_image (j + j0_im[n], i + i0_im[n])
                * (*responses[n])(j + j0_ft[n], i + i0_ft[n]);
            }

        double scale = signal[n] / sum;

        if (iter != 1 && boosting && boost_function)
          {
            scale = boost_function (scale);
          }
        for (int i = 0; i < i1_im[n] - i0_im[n]; ++i)
          for (int j = 0; j < j1_im[n] - j0_im[n]; ++j)
            {
              local_correction (j + j0_im[n], i + i0_im[n])
                += (*responses[n])(j + j0_ft[n], i + i0_ft[n]) * scale;
            }
      }
  }
}

namespace hires
{
void Footprint::compute_correction (
    const std::array<int,2> &nxy, const arma::mat &signal_image,
    const int &iter, const bool &boosting, 
    const std::function<double(double)> &boost_function,
    arma::mat &correction) const
{
  // FIXME: Currently use all available cores.  It would be better to
  // make this configurable.
  const size_t num_threads=std::thread::hardware_concurrency();

  std::vector<arma::mat> local_correction(num_threads);
  std::vector<std::thread> threads;

  for (size_t i = 0; i < num_threads; ++i)
    {
      threads.emplace_back (thread_callback, i, num_threads, signal_image,
                            responses, signal, boosting, boost_function,
                            nxy, iter, j0_im, j1_im, i0_im, i1_im, j0_ft, i0_ft,
                            std::ref (local_correction[i]));
    }
  for (auto &t : threads)
    t.join ();
  
  correction.zeros (nxy[1], nxy[0]);
  for(size_t i=0; i<num_threads; ++i)
    {
      correction+=local_correction[i];
    }
}
}
