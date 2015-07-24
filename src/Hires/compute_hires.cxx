#include <fstream>
#include <boost/math/constants/constants.hpp>
#include <mlpack/methods/lars/lars.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "../Hires.hxx"

void hires::Hires::compute_hires (const boost::filesystem::path &)
// void hires::Hires::compute_hires (const boost::filesystem::path &Drf_file)
{
  // FIXME: hard coded for 44 GHz
  const double sigma_drf2(1.07189706078e-5);
  const double pi(boost::math::constants::pi<double>());

  /// Set up binned data
  double variance;
  arma::vec binned_data;
  size_t bins, min_bins(8), max_bins(64);
  for (size_t b=min_bins; b<=max_bins; b*=2)
    {
      std::vector<boost::accumulators::accumulator_set
                  < double, boost::accumulators::stats
                    < boost::accumulators::tag::count,
                      boost::accumulators::tag::mean,
                      boost::accumulators::tag::variance > > >
        accumulators(b*b);
      /// Offsets are b/2, not (b-1)/2, since the floor function needs
      /// an offset of 0.5.
      const double pix_per_radian=b/(nxy[0]*radians_per_pix),
        offset(b/2.0);
      for (auto &sample: samples)
        {
          double xi ((sample.x * pix_per_radian) + offset);
          double yi ((sample.y * pix_per_radian) + offset);
          if(xi<0 || xi>=b || yi<0 || yi>=b)
            continue;
          int i_int (std::floor(xi)), j_int (std::floor(yi));
          accumulators[i_int + b*j_int](sample.signal);
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
      if (boost::accumulators::count(median) != b*b)
        {
          if (b==min_bins)
            throw hires::Exception
              ("Could not compute the variance of the input.  "
               "The samples must cover every bin in an 8x8 grid");
          else
            break;
        }
      variance=boost::accumulators::median(median);
      bins=b;
      binned_data.resize(b*b);
      for (size_t ii=0; ii<b*b; ++ii)
        binned_data(ii)=boost::accumulators::mean(accumulators[ii]);
    }


  /// Compute transfer function
  const double min_x(0), min_y(0), delta_bin_x(nxy[0]*radians_per_pix/bins),
    delta_bin_y(nxy[1]*radians_per_pix/bins);
  
  arma::mat A(bins*bins,nxy[0]*nxy[1]);
  for(size_t iy=0; iy<nxy[1]; ++iy)
    {
      const double y=min_y + (iy+0.5)*radians_per_pix;
      for(size_t ix=0; ix<nxy[0]; ++ix)
        {
          const double x=min_x + (ix+0.5)*radians_per_pix;

          for(size_t iy_bin=0; iy_bin<bins; ++iy_bin)
            {
              const double bin_y=min_y + (iy_bin+0.5)*delta_bin_y;
              for(size_t ix_bin=0; ix_bin<bins; ++ix_bin)
                {
                  const double bin_x=min_x + (ix_bin+0.5)*delta_bin_x;
                  double dx(x - bin_x), dy(y-bin_y);
                  A(ix_bin + bins*iy_bin,ix + nxy[0]*iy)=
                    radians_per_pix*radians_per_pix
                    *exp(-(dx*dx+dy*dy)/(2*sigma_drf2))/(2*pi*sigma_drf2);
                }
            }
        }
    }
      
  /// LARS
  const double data_scale=std::max(arma::max(binned_data),
                                   std::abs(arma::min(binned_data)));
  std::cout << "variance: " << variance << " "
            << std::sqrt(variance) << " "
            << data_scale << "\n";

  mlpack::regression::LARS lars(true,variance/data_scale,
                                variance/(data_scale*data_scale));
  arma::vec lars_image;
  lars.Regress(A,binned_data,lars_image,false);

  {
    std::ofstream outfile("elastic_net");
    for (size_t x=0; x<nxy[0]; ++x)
      for (size_t y=0; y<nxy[1]; ++y)
        outfile << x << " "
                << y << " "
                << lars_image(x+nxy[0]*y) << " "
                << "\n";
  }

  /// MCM with binning
  const size_t max_iterations=32;
  arma::vec f(nxy[0]*nxy[1]), f_new;
  const double total_flux=arma::mean(binned_data);
  f.fill(total_flux);
  for (size_t iteration=0; iteration<max_iterations; ++iteration)
    {
      arma::vec Af=A*f;
      arma::vec gAf=binned_data/Af;
      f_new=f%(A.t()*gAf);
      f=f_new*(total_flux/arma::mean(f_new));
      {
        std::ofstream outfile("mcm_" + std::to_string(iteration+1));
        for (size_t x=0; x<nxy[0]; ++x)
          for (size_t y=0; y<nxy[1]; ++y)
            outfile << x << " "
                    << y << " "
                    << f(x+nxy[0]*y) << " "
                    << "\n";
      }
    }
}

