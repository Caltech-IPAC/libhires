#include "../Hires.hxx"

void hires::Hires::dump_params ()
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
