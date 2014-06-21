#include "../Hires_Parameters.hxx"

namespace hires {

void Hires_Parameters::set_default_values()
{
    data_type = Data_Type::planck;
    hires_mode = Hires_Mode::both;
    drf_prefix = "";
    ctype1 = "RA -- TAN";
    ctype2 = "DEC -- TAN";
    boost_type = "TIMES_2";
    starting_image = "flat";
    beam_starting_image = "flat";
    flux_units = "??";
    log_filename = "logfile.log";
    outfile_types = {"flux"},
    ni = 500;
    nj = 500;
    boost_max_iter = 0;
    footprints_per_pix = 1;
    beam_spike_n = 5;
    radians_per_pix = 60*boost::math::constants::pi<double>()/(3600*180);
    crval1 = 0.0;
    crval2 = 0.0;
    min_sample_flux = std::numeric_limits<double>::min();
    angle_tolerance = 2.5;
    beam_spike_height = 10;

    std::string kwd("AUTHOR");
    std::string val("LIBHIRES");
    std::string comment("");
    fits_keywords.push_back(std::make_tuple(kwd, val, comment));

    boost_func = [](const double &x) {return x+x-1.0;};

    iter_max = 0;
    iter_list.push_back(0);
}
}
