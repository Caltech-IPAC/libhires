#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{

void Hires::write_file (const std::string &output_prefix,
                        const Image_Type &type)
{
  std::map<Image_Type,std::tuple<std::string,std::string,bool,arma::mat *> >
    image_mapping=
    {{Image_Type::hires_image,
      std::make_tuple("hires", "HIRES image", true, &flux_images)},
     {Image_Type::hires_covariance,
      std::make_tuple("hires_cov", "HIRES covariance", false, &wgt_image)},
     {Image_Type::hires_correction,
      std::make_tuple("hires_cfv", "HIRES correction factor variance", false,
                      &cfv_images)},
     {Image_Type::hires_beam, std::make_tuple("hires_beam", "HIRES beam", false,
                                              &beam_images)},
     {Image_Type::minimap_image, std::make_tuple("minimap", "Minimap Image",
                                                 true, &minimap)},
     {Image_Type::minimap_hitmap,
      std::make_tuple("minimap_hitmap", "Minimap hitmap", false, &hitmap)}};

  std::string filename=output_prefix
    + std::get<0>(image_mapping[type]);
  if (iteration != 0)
    filename+= "_" + std::to_string(iteration);
  filename+= ".fits";

  std::vector<std::pair<std::string, std::pair<std::string, std::string> > >
  file_specific_keywords;

  file_specific_keywords.push_back ({"FILETYPE",
        {std::get<1>(image_mapping[type]), ""}});

  if (iteration!=0)
    {
      file_specific_keywords.push_back
        ({"ITERNUM", {std::to_string(iteration), "HIRES iteration number"}});
              

      file_specific_keywords.push_back
        ({"DRF_IN", {drf_file.string (), "Detector Response File"}});
    }

  /// Is it a flux?
  if (!std::get<2>(image_mapping[type]))
    {
      std::vector<std::pair<std::string, std::pair<std::string, std::string> > >
        sanitized_keywords;
      for(auto &k: file_specific_keywords)
        {
          if(k.first!="BUNIT")
            sanitized_keywords.push_back(k);
        }
      write_fits (*std::get<3>(image_mapping[type]), sanitized_keywords,
                  filename);
    }
  else
    {
      write_fits (*std::get<3>(image_mapping[type]), file_specific_keywords,
                  filename);
    }
}
} // end namespace hires
