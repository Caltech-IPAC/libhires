#include <CCfits>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include <valarray>

#include "../Sample.hxx"
#include "../Gnomonic.hxx"

namespace hires
{
  void read_one_IN_planck(const boost::filesystem::path &filename,
                          const Gnomonic &projection,
                          std::vector<Sample> &samples)
  {
    std::string stem(filename.filename().stem().string());
    auto end(stem.rfind("-"));
    auto begin(stem.rfind("-",end));
    int detector_id(boost::lexical_cast<int>(stem.substr(begin+1,end)));

    std::map<int,std::vector<double> > det_offset_scale_rms
    {{1,{33.858, 10351.1, 1.086}},
        {2,{34.189,  9590.0, 1.097}},
          {3,{33.840, 10640.6, 1.085}}};
    
    std::vector<double> &offset_scale_rms(det_offset_scale_rms[detector_id]);
    double offset(offset_scale_rms[0]);
    double scale(offset_scale_rms[1]);

    CCfits::FITS fits_file(filename.string(),CCfits::Read,true);
    CCfits::ExtHDU& img = fits_file.extension(1); 

    std::vector<std::valarray<double> > phi_vector, theta_vector, signal_vector;
    img.column("PHI").readArrays(phi_vector,0,img.column("PHI").rows());
    img.column("THETA").readArrays(theta_vector,0,img.column("THETA").rows());
    img.column("SIGNAL").readArrays(signal_vector,0,img.column("SIGNAL").rows());

    for(size_t i=0;i<phi_vector.size();++i)
      {
        std::valarray<double> &glon(phi_vector[i]), &theta(theta_vector[i]),
          &signal(signal_vector[i]);
        const size_t n(glon.size());
        std::valarray<double> glat(n), offset_signal(n), x(n), y(n);

        offset_signal=(signal - offset)/scale;
        glat=boost::math::constants::pi<double>()/2 - theta;
        std::tie(x,y)=projection.lonlat2xy(glon,glat);
        // FIXME: Why is this hard coded to 1?
        // samples.emplace_back(-x,y,offset_signal,detector_id);
        samples.emplace_back(-x,y,offset_signal,1);
      }
  }
}

