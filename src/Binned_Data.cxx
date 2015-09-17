#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "Binned_Data.hxx"
#include "Exception.hxx"

hires::Binned_Data::Binned_Data (const std::vector<Sample> &samples,
                                 const std::array<size_t,2> &nxy,
                                 const double &radians_per_pix,
                                 const size_t &max_bins)
{
  boost::accumulators::accumulator_set
    < double, boost::accumulators::stats
      < boost::accumulators::tag::variance > >
    variance_accumulator;
  for (auto &sample: samples)
    variance_accumulator (sample.signal);
  variance=boost::accumulators::variance (variance_accumulator);
    
  size_t min_bins(8);
  /// Keep subdividing until we get empty cells or we get to the max
  /// number of bins.
  for (size_t b=min_bins; b<=max_bins; b*=2)
    {
      std::cout << "bin: " << max_bins << " " << b << "\n";
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
      // variance=boost::accumulators::median(median);
      num_bins=b;
      data.set_size(b*b);
      for (size_t ii=0; ii<accumulators.size(); ++ii)
        data(ii)=boost::accumulators::mean(accumulators[ii]);
    }
}

