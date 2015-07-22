#include "../Hires.hxx"
#include "../Detector.hxx"


#include <random>
#include <fstream>
#include <boost/math/constants/constants.hpp>

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
  const size_t M(51);
  const double delta_x=(max_x-min_x)/(M-1),
    pi(boost::math::constants::pi<double>());

  Eigen::VectorXd d(N), x(N);
  for (size_t i=0; i<N; ++i)
    {
      x(i)=uniform(gen);
      d(i)=exp(-x(i)*x(i)/(2*sigma_detector2))/sqrt(2*pi*sigma_detector2)
        + normal(gen)*noise_scale;
    }
  // std::vector<size_t> bins(101);
  // for (size_t i=0; i<1000000; ++i)
  //   {
  //     size_t index=std::min<size_t>(std::max<size_t>((normal(gen)+10)*100/20,0),100);
  //     ++bins[index];
  //   }
  // for (auto &bin: bins)
  //   std::cout << bin << "\n";
  
  Eigen::MatrixXd A(N,M);
  for(size_t i=0; i<N; ++i)
    for(size_t j=0; j<M; ++j)
      {
        double x_i=min_x + j*delta_x;
        double dx=x(i) - x_i;
        A(i,j)=delta_x*exp(-dx*dx/2)/sqrt(2*pi*sigma_drf2);
      }

  Eigen::MatrixXd AtA=A.transpose()*A;
  Eigen::VectorXd Ag=A.transpose()*d;
  Eigen::VectorXd f;
  
  // Diffusion
  // if(0)
  //   {
  //     Eigen::MatrixXd nabla(M,M);
  //     nabla.setZero();
  //     for(size_t j=1; j<M-1; ++j)
  //       {
  //         nabla(j,j-1)=1;
  //         nabla(j,j)=-2;
  //         nabla(j,j+1)=1;
  //       }
  //     nabla(0,0)=-1;
  //     nabla(0,1)=1;
  //     nabla(M-1,M-2)=1;
  //     nabla(M-1,M-1)=-1;
  //     double f_norm=1; // FIXME: this should compute rms of f or g
  //     double lambda=-1*noise_rms*noise_rms/(f_norm);
  
  //     // std::cout << "nabla\n" << nabla << "\n";
  //     Eigen::MatrixXd AtA=(A.transpose()*A + lambda * nabla);
  //   }
  
  // Tikhonov
  {
    Eigen::MatrixXd nabla(M,M);
    nabla.setZero();
    for(size_t j=0; j<M; ++j)
      nabla(j,j)=1;
    double f_scale=d.mean();
    double lambda=noise_rms*noise_rms*N/(M*f_scale*f_scale);
    AtA+=lambda * nabla;
  
    Eigen::VectorXd rhs=Ag.array()+lambda*f_scale;
    Eigen::VectorXd tikhonov=AtA.ldlt().solve(rhs);
    {
      std::ofstream outfile("tikhonov");
      outfile << tikhonov << "\n";
    }
    std::cout << "tikhonov " << f_scale << " "
              << tikhonov.mean() << "\n";
    f=tikhonov;
  }
  // // Lasso
  // Eigen::MatrixXd nabla(M,M);
  // nabla.setZero();
  // for(size_t j=0; j<M; ++j)
  //   nabla(j,j)=1;
  // double f_scale=1.0;
  // // double f_scale=d.mean();
  // double lambda=noise_rms*noise_rms*0.01;
  // Eigen::VectorXd f=(AtA + lambda*nabla).ldlt().solve(Ag);
  // {
  //   std::ofstream outfile("normal_lsq_0");
  //   outfile << f << "\n";
  // }
  // Eigen::VectorXd tanh_f(M);
  // for (int iteration=1; iteration<10; ++iteration)
  //   {
  //     for(size_t j=0; j<M; ++j)
  //       {
  //         double t_f=std::tanh(f(j)/f_scale);
  //         tanh_f(j)=t_f;
  //         nabla(j,j)=1-t_f*t_f;
  //       }
  //     // FIXME: Normalize f so that we never need to divide by the
  //     // noise?
  //     Eigen::VectorXd df=(AtA + lambda*nabla).ldlt().solve(Ag - AtA*f - lambda*tanh_f);
  //     f+=df;
  //     // if(iteration%10==0)
  //     {
  //       std::ofstream outfile("normal_lsq_" + std::to_string(iteration));
  //       outfile << f << "\n";
  //     }


  //     for(size_t j=0; j<M; ++j)
  //       {
  //         double t_f=std::tanh(f(j)/f_scale);
  //         tanh_f(j)=t_f;
  //       }
  //     std::cout << "||Af-d||^2 " << iteration << " "
  //               << (Ag - AtA*f - lambda*tanh_f).norm() << " "
  //               << (Ag - AtA*f).norm() << " "
  //               << noise_rms
  //               << "\n";
  //   }

  // // Richardson-Lucy
  // // This is so bogus
  // Eigen::VectorXd lucy(M);
  // lucy.setConstant(d.mean());
  // const double K=lucy.sum()/d.sum();
  // {
  //   std::ofstream outfile("lucy_0");
  //   outfile << lucy << "\n";
  // }
  // std::cout << "||Af-d||^2 " << "0" << " " << sqrt((A*lucy-d).squaredNorm()) << " "
  //           << noise_rms
  //           << "\n";

  // for (int iteration=1; iteration<10; ++iteration)
  //   {
  //     Eigen::VectorXd Af=A*lucy;
  //     Eigen::VectorXd resid=d.cwiseQuotient(Af);
  //     Eigen::VectorXd factor=A.transpose()*resid;
  //     // std::cout << "lucy\n"
  //     //           << lucy
  //     //           << "\nK\n" << K
  //     //           << "\nresid\n" << resid
  //     //           << "\nAf\n" << Af
  //     //           // << "\nx\n" << x
  //     //           << "\nd\n" << d
  //     //           << "\nfactor\n" << factor << "\n";
  //     lucy=lucy.cwiseProduct(K*factor);
  //     if(iteration%100==0)
  //       {
  //         std::ofstream outfile("lucy_" + std::to_string(iteration/100));
  //         outfile << lucy << "\n";
  //       }
  //     std::cout << "||Af-d||^2 " << iteration << " " << (A*lucy-d).norm() << " "
  //               << noise_rms
  //               << "\n";
  //   }

  // Eigen::VectorXd f=AtA.ldlt().solve(Ag);
  // {
  //   std::ofstream outfile("normal_lsq");
  //   outfile << f << "\n";
  // }

  // // Eigen::VectorXd f;
  // f=A.householderQr().solve(d);
  // {
  //   std::ofstream outfile("HouseholderQr");
  //   outfile << f << "\n";
  // }
  // f=A.colPivHouseholderQr().solve(d);
  // {
  //   std::ofstream outfile("colPivHouseholderQr");
  //   outfile << f << "\n";
  // }
  // {
  //   std::ofstream outfile("fullPivHouseholderQr");
  //   outfile << A.fullPivHouseholderQr().solve(d) << "\n";
  // }
  // {
  //   std::ofstream outfile("svd");
  //   outfile << A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(d) << "\n";
  // }

  {
    std::ofstream outfile("x_dx");
    Eigen::VectorXd Af=A*f;
    for (size_t i=0; i<N; ++i)
      outfile << x(i) << " " << d(i) << " "
              << Af(i) << "\n";
  }
  
  // {
  //   std::ofstream outfile("residual");
  //   outfile << (AtA*f - Ag) << "\n";
  // }

  // {
  //   std::ofstream outfile("AF");
  //   outfile << A*f << "\n";
  // }

  // {
  //   std::ofstream outfile("d");
  //   outfile << d << "\n";
  // }

  // {
  //   std::ofstream outfile("AAF");
  //   outfile << AtA*f << "\n";
  // }

  // Eigen::VectorXd jacobian(M);
  // jacobian.setZero();
  // for (size_t j=0; j<M; ++j)
  //   for (size_t i=0; i<N; ++i)
  //     jacobian(j)+=A(i,j)*A(i,j);

  // f.resize(M);
  // f.setZero();
  // for (size_t iteration=0; iteration<10; ++iteration)
  //   {
  //     Eigen::VectorXd Afd=A*f - d;
  //     Eigen::VectorXd resid=A.transpose()*Afd;
  //     f-=resid.cwiseQuotient(jacobian);
  //     // if(iteration%100==0)
  //       std::cout << "iteration: " << iteration << "\n"
  //                 << f << "\n";
  //   }


  // for (size_t i=0; i<f.size(); ++i)
  //   std::cout << f[i] << "\n";
  
  // for (auto &sample: samples)
  //   {
  //     for (size_t y=0; y<nxy[1]; ++y)
  //       for (size_t x=0; x<nxy[0]; ++x)
  //         {
  //           double xi=(x-i_offset)*radians_per_pix;
  //           double yi=(y-j_offset)*radians_per_pix;
            
  //           // FIXME: This ignores rotation.  Some of that is handled
  //           // in Footprint, but I am still unsure what footprint
  //           // realy does.
  //           double response=d.response(sample.x - xi, sample.y - yi);
  //           Ag(x + nxy[0]*y)+=response * sample.signal;

  //           for (size_t yt=0; yt<nxy[1]; ++yt)
  //             for (size_t xt=0; xt<nxy[0]; ++xt)
  //               {
  //                 double xti=(xt-i_offset)*radians_per_pix;
  //                 double yti=(yt-j_offset)*radians_per_pix;

  //                 AA(xt+nxy[0]*yt,x+nxy[0]*y)+=response
  //                   * d.response(sample.x - xti, sample.y - yti);
  //               }
  //         }
  //     if (s%1000==0)
  //       std::cout << "computing sample " << s << " of " << samples.size() << " "
  //                 << nxy[0] << " " << nxy[1] << " "
  //                 << "\n";
  //     ++s;
  //   }
  // std::cout << "solving\n";
  // Eigen::VectorXd answer=AA.transpose().ldlt().solve(Ag);
  // std::cout << "solved\n";
  // hires.resize (nxy[0],nxy[1]);

  // for(size_t y=0; y<nxy[1]; ++y)
  //   for(size_t x=0; x<nxy[0]; ++x)
  //     hires(x,y)=answer(x+nxy[0]*y);
}

