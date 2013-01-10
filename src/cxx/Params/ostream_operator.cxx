#include "../Params.hxx"
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>

std::ostream& operator<<(std::ostream& out, const Params &p)
{
  out << "\nInput data file options:"
      << "\n  INFILE_PREFIX " << p.infile_prefix
      << "\n  STARTING_IMAGE " << p.starting_image
      << boost::format("\n  MIN_SAMPLE_FLUX %.6f") % p.min_sample_flux

      << "\n\nDRF (detector response files) to use:"
      << "\n  DRF_PREFIX " << p.drf_prefix

      << "\n\nOutput image geometry:"
      << "\n  NPIX " << p.ni << " " << p.nj
      << boost::format("\n  DEG_PER_PIX %.6f")
    % (p.radians_per_pix*180/boost::math::constants::pi<double>())
      << boost::format("\n  CRVAL1 %.5f") % p.crval1
      << boost::format("\n  CRVAL2 %.5f") % p.crval2
      << "\n  CTYPE1 " << p.ctype1
      << "\n  CTYPE2 " << p.ctype2

      << "\n\nOutput file options:"
      << "\n  OUTFILE_PREFIX " << p.outfile_prefix
      << "\n  OUTFILE_TYPES [";
  for(auto &o: p.outfile_types)
    out << "'" << o << "' ";
  out << "]\n  ITER_MAX " << p.iter_max
      << "\n  ITER_LIST [";
  for(auto &o: p.iter_list)
    out << o << " ";
  out << "]\n  FLUX_UNITS " << p.flux_units;

  if(find(p.outfile_types.begin(),p.outfile_types.end(),"beam")
     !=p.outfile_types.end())
    {
      out << "\n\nBeam image file options:"
          << "\n  BEAM_SPIKE_N " << p.beam_spike_n
          << "\n  BEAM_SPIKE_HEIGHT " << p.beam_spike_height
          << "\n  BEAM_STARTING_IMAGE " << p.beam_starting_image;
    }
  out << "\n\nAdditional FITS keywords:";
  for(auto &k: p.fits_keywords)
    out << "\n  KWD ('" << std::get<0>(k) << "', '" << std::get<1>(k) << "', '" << std::get<2>(k) << "')";

  out << "\n\nAccelerated correction option:";
  if(p.boost_max_iter>0)
    out << "\n  BOOST " << p.boost_type
        << " for ITER 2 to " << p.boost_max_iter;
  else
    out << "\n  BOOST (none)";

  out << "\n\nFootprint accuracy options:"
      << boost::format("\n  ANGLE_TOLERANCE %.2f") % p.angle_tolerance
      << "\n  FOOTPRINTS_PER_PIX " << p.footprints_per_pix
      << "\n";
  return out;
}
