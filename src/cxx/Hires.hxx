#ifndef HIRES_PARAMS_HXX
#define HIRES_PARAMS_HXX

#include <string>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <functional>
#include <armadillo>
#include <map>
#include "Sample.hxx"

namespace hires
{
  class Hires
  {
  public:
    std::string drf_prefix, ctype1, ctype2,
      boost_type, starting_image, beam_starting_image, flux_units, log_filename;
    std::vector<std::string> outfile_types;
    int ni, nj, iter_max, boost_max_iter, footprints_per_pix, beam_spike_n;
    double radians_per_pix, crval1, crval2, min_sample_flux, angle_tolerance,
      beam_spike_height;

    arma::mat hitmap, minimap;

    std::vector<int> iter_list;
    std::vector<std::tuple<std::string,std::string,std::string>> fits_keywords;
    std::function<double (double)> boost_func;

    Hires(const std::string &Data_type, const std::string &Hires_mode, const std::vector<std::string> &param_files);


    enum class Data_Type {planck, spire} data_type;
    enum class Hires_Mode {hires, minimap, both} hires_mode;

    void compute_images(std::vector<Sample> &samples,
                        const std::string &outfile_prefix)
    {
      arma::mat wgt_image;
      std::map<int,arma::mat> flux_images;
      std::map<int,arma::mat> cfv_images;
      std::map<int,arma::mat> beam_images;
      compute_images(wgt_image,flux_images,cfv_images,
                     beam_images,samples,outfile_prefix);
    }
    void compute_images(arma::mat &wgt_image,
                        std::map<int,arma::mat> &flux_images,
                        std::map<int,arma::mat> &cfv_images,
                        std::map<int,arma::mat> &beam_images,
                        std::vector<Sample> &samples)
    {
      compute_images(wgt_image,flux_images,cfv_images,
                     beam_images,samples,"");
    }
    void compute_images(arma::mat &wgt_image,
                        std::map<int,arma::mat> &flux_images,
                        std::map<int,arma::mat> &cfv_images,
                        std::map<int,arma::mat> &beam_images,
                        std::vector<Sample> &samples,
                        const std::string &outfile_prefix);

    void iterate(arma::mat &wgt_image,
                        std::map<int,arma::mat> &flux_images,
                        std::map<int,arma::mat> &cfv_images,
                        std::map<int,arma::mat> &beam_images,
                        std::vector<Sample> &samples,
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
