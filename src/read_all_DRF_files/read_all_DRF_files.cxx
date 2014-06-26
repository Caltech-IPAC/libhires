#include "../Hires.hxx"
#include "../Detector.hxx"

namespace hires
{
  std::map<int,Detector>
  read_all_DRF_planck(const std::string &DRF_prefix);

  std::map<int,Detector>
  read_all_DRF_files(const std::string &DRF_prefix)
  {
      return read_all_DRF_planck(DRF_prefix);
  }
}
