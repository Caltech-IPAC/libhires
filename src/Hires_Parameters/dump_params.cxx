#include "../Hires_Parameters.hxx"

namespace hires
{

void Hires_Parameters::dump_params ()
{
  std::cout << "drf_prefix " << drf_prefix.c_str () << "\n";
  std::cout << "ctype1 " << ctype1.c_str () << "\n";
  std::cout << "ctype2 " << ctype2.c_str () << "\n";
  std::cout << "boost_type " << boost_type.c_str () << "\n";
  std::cout << "starting_image  " << starting_image.c_str () << "\n";
  std::cout << "beam_starting_image " << beam_starting_image.c_str () << "\n";
  std::cout << "flux_units " << flux_units.c_str () << "\n";
  std::cout << "log_filename " << log_filename.c_str () << "\n";
  std::cout << "ni " << ni << "\n";
  std::cout << "nj " << nj << "\n";
  std::cout << "boost_max_iter " << boost_max_iter << "\n";
  std::cout << "footprints_pr_pix " << footprints_per_pix << "\n";
  std::cout << "beam_spike_n " << beam_spike_n << "\n";
  std::cout << "radians_per_pix " << radians_per_pix << "\n";
  std::cout << "crval1 " << crval1 << "\n";
  std::cout << "crval2 " << crval2 << "\n";
  std::cout << "min_sample_flux " << min_sample_flux << "\n";
}
}

/*
        Hires(hp.data_type,
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
           hp.boost_func);
*/
