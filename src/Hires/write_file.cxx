#include "../Hires.hxx"

namespace hires
{

void Hires::write_file (const std::string &output_prefix,
                        const std::string &filename,
                        const std::string &filetype,
                        const arma::mat &image)
{
  std::string file;
  if(output_prefix=="-")
    {
      file="-";
    }
  else
    {
      file=output_prefix + filename + ".fits";
      boost::filesystem::path output_path=boost::filesystem::path(output_prefix+filename);
      if (!exists(output_path.parent_path()))
        create_directories(output_path.parent_path());
    }

  std::vector<std::pair<std::string, std::pair<std::string, std::string> > >
  file_specific_keywords;

  file_specific_keywords.push_back ({"FILETYPE", {filetype, ""}});
  write_fits (image, file_specific_keywords, file);
}
}
