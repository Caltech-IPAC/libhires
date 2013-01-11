#include <boost/filesystem.hpp>
#include <vector>
#include "../Sample.hxx"
#include "../Params.hxx"
#include "../Gnomonic.hxx"
#include "../Exception.hxx"

namespace hires
{
  void read_one_IN_planck(const boost::filesystem::path &p,
                          const Gnomonic &projection,
                          std::vector<Sample> &samples);

  std::vector<Sample> read_all_IN_files(const Params::Data_Type &dt,
                                        const std::string &prefix,
                                        const Gnomonic &projection)
  {
    std::vector<Sample> samples;
    boost::filesystem::path p(prefix);
    boost::filesystem::directory_iterator end;
    for(boost::filesystem::directory_iterator d(p.parent_path());
        d!=end; ++d)
      {
        std::string s(p.filename().string()), s2(d->path().filename().string());

        if((s==s2.substr(0,s.size())) && s2.size()>5
           && s2.substr(s2.size()-5)==".fits")
          {
            if(dt==Params::Data_Type::planck)
              read_one_IN_planck(d->path(),projection,samples);
            else
              abort();
          }
      }
    if(samples.empty())
      throw Exception("No input files for: " + prefix);
    return samples;
  }
}
