#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "../Hires.hxx"

hires::Hires::Binned_Data hires::Hires::bin_data ()
{
  Binned_Data binned_data;
  size_t min_bins(8), max_bins(64);
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
      binned_data.variance=boost::accumulators::median(median);
      binned_data.num_bins=b;
      binned_data.data.set_size(b*b);
      for (size_t ii=0; ii<accumulators.size(); ++ii)
        binned_data.data(ii)=boost::accumulators::mean(accumulators[ii]);
    }
  return binned_data;
}

