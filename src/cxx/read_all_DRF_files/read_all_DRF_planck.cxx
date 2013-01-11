#include "../Detector.hxx"

namespace hires
{
  std::map<int,Detector>
  read_all_DRF_planck(const std::string &DRF_prefix)
  {
    std::map<int,Detector> detectors;

    boost::filesystem::path prefix(DRF_prefix);

    for(boost::filesystem::directory_iterator p(prefix.parent_path());
        p!=boost::filesystem::directory_iterator(); ++p)
      {
        if(p->path().string().substr(0,prefix.string().size())==prefix.string())
          {
            Detector d(p->path());
            /* A trick so that the underlying array is not copied */
            std::swap(detectors[d.id],d);
          }
      }
    if(detectors.empty())
      throw Exception("No valid detectors found in " + DRF_prefix);
    return detectors;
  }
}
