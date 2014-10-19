#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{

void Hires::write_file (arma::mat image, std::string filename,
                        const char *desc, int iter, int isflux)
{
  std::vector<std::pair<std::string, std::pair<std::string, std::string> > >
  file_specific_keywords;

  file_specific_keywords.push_back ({"FILETYPE",{desc, ""}});
        

  if (iter!=0)
    {
      file_specific_keywords.push_back
        ({"ITERNUM", {std::to_string(iter), "HIRES iteration number"}});
              

      file_specific_keywords.push_back
        ({"DRF_IN", {drf_file.string (), "Name of Detector Response Files"}});
    }

  if (!isflux)
    {
      std::vector<std::pair<std::string, std::pair<std::string, std::string> > >
        sanitized_keywords;
      for(auto &k: file_specific_keywords)
        {
          if(k.first!="BUNIT")
            sanitized_keywords.push_back(k);
        }
      write_fits (image, sanitized_keywords, filename);
    }
  else
    {
      write_fits (image, file_specific_keywords, filename);
    }
}
} // end namespace hires
