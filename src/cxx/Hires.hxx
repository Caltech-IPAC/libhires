#ifndef HIRES_PARAMS_HXX
#define HIRES_PARAMS_HXX

#include <string>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <functional>
#include <armadillo>
#include <map>
#include "Hires_Parameters.hxx"
#include "Sample.hxx"
#include "Exception.hxx"

namespace hires {

class Hires : public Hires_Parameters
  {
  public:
    arma::mat hitmap, minimap;


// Constructors
// takes Hires_Parameters as single argument
    Hires(const Hires_Parameters &hp) : Hires_Parameters(hp.data_type, 
           hp.hires_mode, 
           hp.drf_prefix,
           hp.ctype1,
           hp.ctype2,
           hp.boost_type,
           hp.starting_image,
           hp.beam_starting_image,
           hp.flux_units,
           hp.log_filename,
           hp.outfile_types,
           hp.ni,
           hp.nj,
           hp.boost_max_iter,
           hp.footprints_per_pix,
           hp.beam_spike_n,
           hp.radians_per_pix,
           hp.crval1,
           hp.crval2,
           hp.min_sample_flux,
           hp.angle_tolerance,
           hp.beam_spike_height,  
           hp.fits_keywords,
           hp.boost_func,
           hp.iter_max,
           hp.iter_list) {
    };

    void iterate(arma::mat &wgt_image,
                        std::map<int,arma::mat> &flux_images,
                        std::map<int,arma::mat> &cfv_images,
                        std::map<int,arma::mat> &beam_images,
                        std::vector<hires::Sample> &samples,
                        int &iter);

    void write_output(arma::mat &wgt_image,
                        std::map<int,arma::mat> &flux_images,
                        std::map<int,arma::mat> &cfv_images,
                        std::map<int,arma::mat> &beam_images,
                        int &iter,
                        const std::string &outfile_prefix);


    arma::mat spike_image();
    arma::mat start_image(const std::string &filename, int &iter_start);

    void write_fits(const arma::mat &image, const std::string &file_type,
                    const int &iter, const std::string &outfile_prefix);

    void write_fits(const arma::mat &image, const std::string &file_type,
                    const std::string &outfile_prefix)
    {
      write_fits(image,file_type,std::numeric_limits<int>::max(),outfile_prefix);
    }

    static std::string type_string(const Data_Type &dt)
    {
      if(dt==Data_Type::planck)
        return "planck";
      else if(dt==Data_Type::spire)
        return "spire";
      abort();
      return "";
    }
  };

  std::ostream& operator<<(std::ostream& out, const Hires &p);
}
#endif
