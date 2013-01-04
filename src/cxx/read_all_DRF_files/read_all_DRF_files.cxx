#include "../Params.hxx"
#include "../Detector.hxx"
#include "../logger.hxx"

std::map<std::string,Detector>
read_all_DRF_planck(const std::string &DRF_prefix);

std::map<std::string,Detector>
read_all_DRF_files(const Params::Data_Type &dt, const std::string &DRF_prefix)
{
  switch(dt)
    {
    case Params::Data_Type::planck:
      return read_all_DRF_planck(DRF_prefix);
      break;
    // case Params::Data_Type::spire:
    //   return read_all_DRF_spire(DRF_prefix,logger);
    //   break;
    default:
      LOG4CXX_FATAL(logger,"No valid DRF reader known for " << Params::type_string(dt) << "\n");
      break;
    }
  return std::map<std::string,Detector>();
}
