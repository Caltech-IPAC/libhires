#include "../Detector.hxx"

namespace hires
{
std::map<int, Detector> read_all_DRF_planck (const boost::filesystem::path &DRF_prefix)
{
  std::map<int, Detector> detectors;

  Detector d (DRF_prefix);
  /* A trick so that the underlying array is not copied */
  std::swap (detectors[d.id], d);


  // for (boost::filesystem::directory_iterator p (DRF_prefix.parent_path ());
  //      p != boost::filesystem::directory_iterator (); ++p)
  //   {
  //     if (p->path ().string ().substr (0, DRF_prefix.string ().size ())
  //         == DRF_prefix.string ())
  //       {
  //         Detector d (p->path ());
  //         /* A trick so that the underlying array is not copied */
  //         std::swap (detectors[d.id], d);
  //       }
  //   }
  // if (detectors.empty ())
  //   throw Exception ("No valid detectors found in " + DRF_prefix.string());
  return detectors;
}
}
