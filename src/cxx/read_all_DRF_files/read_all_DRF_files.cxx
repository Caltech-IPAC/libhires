#include "../Hires.hxx"
#include "../Detector.hxx"

namespace hires
{
  std::map<int,Detector>
  read_all_DRF_planck(const std::string &DRF_prefix);

  std::map<int,Detector>
  read_all_DRF_files(const Hires::Data_Type &dt, const std::string &DRF_prefix)
  {
    switch(dt)
      {
      case Hires::Data_Type::planck:
        return read_all_DRF_planck(DRF_prefix);
        break;
        // case Hires::Data_Type::spire:
        //   return read_all_DRF_spire(DRF_prefix);
        //   break;
      default:
        throw Exception("No valid DRF reader known for "
                        + Hires::type_string(dt));
        break;
      }
    return std::map<int,Detector>();
  }
}
