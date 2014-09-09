#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{

void Hires::write_file (arma::mat image, std::string filename,
                        const char *desc, int iter, int isflux)
{
  std::vector<std::tuple<std::string, std::string, std::string> >
  file_specific_keywords;

  file_specific_keywords.push_back (std::make_tuple ("FILETYPE", desc, ""));

  if (iter!=0)
    {
      file_specific_keywords.push_back
        (std::make_tuple ("ITERNUM", std::to_string(iter),
                          "HIRES iteration number"));

      file_specific_keywords.push_back
        (std::make_tuple ("DRF_IN",
                          boost::filesystem::path (drf_prefix).filename ().string ()
                          + "*",
                          "Name of Detector Response Files"));
    }

  if (isflux)
    file_specific_keywords.push_back (std::make_tuple ("BUNIT", flux_units, ""));

  write_fits (image, file_specific_keywords, filename);
}
} // end namespace hires
