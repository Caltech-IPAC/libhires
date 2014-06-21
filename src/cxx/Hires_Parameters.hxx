#ifndef HIRES_PARAMETERS_HXX
#define HIRES_PARAMETERS_HXX

#include <string>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <functional>
#include <armadillo>
#include <map>
#include <boost/math/constants/constants.hpp>

namespace hires {

class Hires_Parameters {

public:
// Enums
    enum class Data_Type {planck, spire} data_type;
    enum class Hires_Mode {hires, minimap, both} hires_mode;

// Elements
    std::string drf_prefix, ctype1, ctype2,
      boost_type, starting_image, beam_starting_image, flux_units, log_filename;
    std::vector<std::string> outfile_types;
    int ni, nj, boost_max_iter, footprints_per_pix, beam_spike_n;
    double radians_per_pix, crval1, crval2, min_sample_flux, angle_tolerance,
      beam_spike_height;

    std::vector<std::tuple<std::string,std::string,std::string>> fits_keywords;
    std::function<double (double)> boost_func;

// Iteration specifics only used by driver, not within Hires_Parameters or Hires
    int iter_max;
    std::vector<int> iter_list;

// Methods

// Constructor #1: Specifies each parameter
    Hires_Parameters(Data_Type data_type, 
        Hires_Mode hires_mode, 
        std::string drf_prefix, 
        std::string ctype1, 
        std::string ctype2,
        std::string boost_type, 
        std::string starting_image,
        std::string beam_starting_image,
        std::string flux_units, 
        std::string log_filename,
        std::vector<std::string> outfile_types,
        int ni, 
        int nj, 
        int boost_max_iter, 
        int footprints_per_pix, 
        int beam_spike_n,
        double radians_per_pix, 
        double crval1, 
        double crval2, 
        double min_sample_flux,
        double angle_tolerance,
        double beam_spike_height,
        std::vector<std::tuple<std::string,std::string,std::string>> fits_keywords,
        std::function<double (double)> boost_func,
        int iter_max,
        std::vector<int> iter_list) :
// Populate all elements
        data_type(data_type),
        hires_mode(hires_mode),
        drf_prefix(drf_prefix),
        ctype1(ctype1),
        ctype2(ctype2),
        boost_type(boost_type),
        starting_image(starting_image),
        beam_starting_image(beam_starting_image),
        flux_units(flux_units),
        log_filename(log_filename),
        outfile_types(outfile_types),
        ni(ni),
        nj(nj),
        boost_max_iter(boost_max_iter),
        footprints_per_pix(footprints_per_pix),
        beam_spike_n(beam_spike_n),
        radians_per_pix(radians_per_pix),
        crval1(crval1),
        crval2(crval2),
        min_sample_flux(min_sample_flux),
        angle_tolerance(angle_tolerance),
        beam_spike_height(beam_spike_height),
        fits_keywords(fits_keywords),
        boost_func(boost_func),
        iter_max(iter_max), 
        iter_list(iter_list) { 
        }

// Constructor #2: Overrides default values based on command line arguments
    Hires_Parameters(const std::string &Data_type, const std::string &Hires_mode,
       const std::vector<std::string> param_str) {
        set_default_values();
        parse_command_line(Data_type, Hires_mode, param_str);
    }

// Populates parameter values based on defaults
    void set_default_values();

// Populates parameter values based on command line arguments
    void parse_command_line(const std::string &Data_type, const std::string &Hires_mode,
       const std::vector<std::string> param_str);

// Diagnostic output
    void dump_params();
};
}
#endif

