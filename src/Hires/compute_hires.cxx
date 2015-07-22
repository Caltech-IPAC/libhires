#include "../Hires.hxx"
#include "../Detector.hxx"


#include <random>
#include <fstream>
#include <boost/math/constants/constants.hpp>
#include <mlpack/methods/lars/lars.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>

void hires::Hires::compute_hires (const boost::filesystem::path &)
// void hires::Hires::compute_hires (const boost::filesystem::path &Drf_file)
{
  // drf_file=Drf_file;
  // Detector d(drf_file);

  // Eigen::MatrixXd AA(nxy[0]*nxy[1],nxy[0]*nxy[1]);
  // AA.setZero ();

  // Eigen::VectorXd Ag(nxy[0]*nxy[1]);
  // Ag.setZero ();

  // double i_offset (nxy[0] / 2.0), j_offset (nxy[1] / 2.0);
  // std::size_t s=0;

  /// Compute the variance by looking at the median variance with
  /// successively finer and finer bins.  Choose the variance when any
  /// of the bins have 0 or 1 sample in it.

  /// FIXME: This is incredibly hokey, but seems to give reasonable
  /// results.  We should really use a more sophisticated statistical
  /// method.
  
  double variance;
  size_t min_bins(8), max_bins(8192);
  for (size_t bins=min_bins; bins<=max_bins; bins*=2)
    {
      std::vector<boost::accumulators::accumulator_set
                  < double, boost::accumulators::stats
                    < boost::accumulators::tag::count,
                      boost::accumulators::tag::variance > > >
        accumulators(bins*bins);
      const double pix_per_radian=bins/(nxy[0]*radians_per_pix),
        offset(bins/2.0);
      for (auto &sample: samples)
        {
          double xi ((sample.x * pix_per_radian) + offset);
          double yi ((sample.y * pix_per_radian) + offset);
          if(xi<0 || xi>=bins || yi<0 || yi>=bins)
            continue;
          int i_int (std::floor(xi)), j_int (std::floor(yi));
          accumulators[i_int + bins*j_int](sample.signal);
        }
      boost::accumulators::accumulator_set
        < double, boost::accumulators::stats
          < boost::accumulators::tag::count,
            boost::accumulators::tag::median> > median;
      for(auto &accumulator: accumulators)
        {
          if(boost::accumulators::count(accumulator)>1)
            median(boost::accumulators::variance(accumulator));
        }
      if (boost::accumulators::count(median) != bins*bins)
        {
          if (bins==min_bins)
            throw hires::Exception
              ("Could not compute the variance of the input.  "
               "The samples must cover every bin in an 8x8 grid");
          else
            break;
        }
      variance=boost::accumulators::median(median);
    }

  std::cout << "variance: " << variance << "\n";

  

  
  std::random_device rd;
  std::mt19937 gen(rd());
  const double min_x(-5), max_x(5);
  std::uniform_real_distribution<> uniform(min_x, max_x);
  const double noise_rms=0.001, noise_scale=1.0;
  std::normal_distribution<> normal(0,noise_rms);
  const double sigma_drf2(1), sigma_signal2(1.0/16),
    sigma_detector2(sigma_drf2 + sigma_signal2);
  
  const size_t N(1000);
  const size_t M(500);
  const double delta_x=(max_x-min_x)/M,
    pi(boost::math::constants::pi<double>());

  arma::vec d(N), x(N);
  for (size_t i=0; i<N; ++i)
    {
      x(i)=uniform(gen);
      d(i)=exp(-x(i)*x(i)/(2*sigma_detector2))/sqrt(2*pi*sigma_detector2)
        + normal(gen)*noise_scale;
    }
  
  arma::mat A(N,M);
  arma::vec x_f(M);
  for(size_t j=0; j<M; ++j)
      {
        x_f(j)=min_x + (j+0.5)*delta_x;
        for(size_t i=0; i<N; ++i)
          {
            double dx=x(i) - x_f(j);
            A(i,j)=delta_x*exp(-dx*dx/(2*sigma_drf2))/sqrt(2*pi*sigma_drf2);
          }
      }

  // Lars
  const double f_scale=std::max(arma::max(d),std::abs(arma::min(d)));
  mlpack::regression::LARS lars(true,noise_rms/f_scale,
            noise_rms*noise_rms/(f_scale*f_scale));
  arma::vec f;
  lars.Regress(A,d,f,false);

  std::cout << "f_scale: " << f_scale << "\n";
  {
    std::ofstream outfile("lars");
    for (size_t m=0; m<M; ++m)
      outfile << x_f(m) << " "
              << f(m) << "\n";
  }

  {
    std::ofstream outfile("x_dx");
    arma::vec Af=A*f;
    for (size_t i=0; i<N; ++i)
      outfile << x(i) << " " << d(i) << " "
              << Af(i) << "\n";
  }
}

