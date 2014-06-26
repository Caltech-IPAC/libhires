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

enum class Image_Type {all, hires, cov, cfv, beam, minimap, hitmap};

class Hires_Parameters {

public:
// Elements
    std::string drf_prefix, ctype1, ctype2,
      boost_type, starting_image, beam_starting_image, flux_units, log_filename;
    std::vector<std::string> outfile_types;
    int ni, nj, boost_max_iter, footprints_per_pix, beam_spike_n;
    double radians_per_pix, crval1, crval2, min_sample_flux, angle_tolerance,
      beam_spike_height;

    std::vector<std::tuple<std::string,std::string,std::string>> fits_keywords;
    std::function<double (double)> boost_func;

// Methods

// Constructor #1: Specifies each parameter
    Hires_Parameters(
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
        std::function<double (double)> boost_func) :
// Populate all elements
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
        boost_func(boost_func) {
        }

// Constructor #2: Overrides default values based on command line arguments
    Hires_Parameters(const std::vector<std::string> param_str) {
        set_default_values();
        parse_command_line(param_str);
    }

// Populates parameter values based on defaults
    void set_default_values();

// Populates parameter values based on command line arguments
    void parse_command_line(const std::vector<std::string> param_str);

// Diagnostic output
    void dump_params();
};
}
#endif

