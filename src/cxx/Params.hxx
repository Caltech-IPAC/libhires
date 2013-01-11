#ifndef HIRES_PARAMS_HXX
#define HIRES_PARAMS_HXX

#include <string>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <functional>
#include <armadillo>
#include <map>

namespace hires
{
  class Params
  {
  public:
    std::string infile_prefix, outfile_prefix, drf_prefix, ctype1, ctype2,
      boost_type, starting_image, beam_starting_image, flux_units, log_filename;
    std::vector<std::string> outfile_types;
    int ni, nj, iter_max, boost_max_iter, footprints_per_pix, beam_spike_n;
    double radians_per_pix, crval1, crval2, min_sample_flux, angle_tolerance,
      beam_spike_height;
    std::vector<int> iter_list;
    std::vector<std::tuple<std::string,std::string,std::string>> fits_keywords;
    std::function<double (double)> boost_func;

    Params(const std::string &Data_type, const std::string &Infile_prefix,
           const std::string &Outfile_prefix,
           const std::vector<std::string> &param_files);

    enum class Data_Type {planck, spire} data_type;

    void compute_images()
    {
      arma::mat wgt_image;
      std::map<int,arma::mat> flux_images;
      std::map<int,arma::mat> cfv_images;
      std::map<int,arma::mat> beam_images;
      compute_images(wgt_image,flux_images,cfv_images,beam_images,true);
    }
    void compute_images(arma::mat &wgt_image,
                        std::map<int,arma::mat> &flux_images,
                        std::map<int,arma::mat> &cfv_images,
                        std::map<int,arma::mat> &beam_images)
    {
      compute_images(wgt_image,flux_images,cfv_images,beam_images,false);
    }
    void compute_images(arma::mat &wgt_image,
                        std::map<int,arma::mat> &flux_images,
                        std::map<int,arma::mat> &cfv_images,
                        std::map<int,arma::mat> &beam_images,
                        const bool &write_images);

    arma::mat spike_image();
    arma::mat start_image(const std::string &filename, int &iter_start);

    void write_fits(const arma::mat &image, const std::string &file_type,
                    const int &iter);

    void write_fits(const arma::mat &image, const std::string &file_type)
    {
      write_fits(image,file_type,std::numeric_limits<int>::max());
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

  std::ostream& operator<<(std::ostream& out, const Params &p);
}
#endif
