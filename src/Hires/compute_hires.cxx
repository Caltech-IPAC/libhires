#include "../Hires.hxx"
#include "../Detector.hxx"


#include <random>
#include <fstream>
#include <boost/math/constants/constants.hpp>
#include <mlpack/methods/lars/lars.hpp>

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

