#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{

void Hires::write_file (const std::string &output_prefix,
                        const std::string &filename,
                        const std::string &filetype,
                        const bool add_drf_filename,
                        const Eigen::MatrixXd &image)
{
  std::string file;
  if(output_prefix=="-")
    {
      file="-";
    }
  else
    {
      file=output_prefix + filename + ".fits";
    }

  std::vector<std::pair<std::string, std::pair<std::string, std::string> > >
  file_specific_keywords;

  file_specific_keywords.push_back ({"FILETYPE", {filetype, ""}});

  if (add_drf_filename)
    file_specific_keywords.push_back
      ({"DRF_IN", {drf_file.string (), "Detector Response File"}});

  write_fits (image, file_specific_keywords, file);
}
}
