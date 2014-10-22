#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{

void Hires::write_output (const std::string &output_prefix)
{
  std::vector<Image_Type> types;
  if(iteration == 0)
    types={Image_Type::hires_covariance, Image_Type::minimap_image,
           Image_Type::minimap_hitmap};
  else
    types={Image_Type::hires_image};

  for(auto &type: types)
    if (output_types.find(type)!=output_types.end())
      write_file (output_prefix, type);
}
}
