#include "../Hires.hxx"

namespace hires
{

void Hires::write_file (const boost::filesystem::path &output_prefix,
                        const std::string &filename,
                        const std::string &filetype,
                        const arma::mat &image)
{
  boost::filesystem::path file(output_prefix);
  if(file!="-")
    {
      file+=filename + ".fits";
      if (!exists(file.parent_path()) && !file.parent_path().empty())
        create_directories(file.parent_path());
    }

  std::vector<std::pair<std::string, std::pair<std::string, std::string> > >
  file_specific_keywords;

  file_specific_keywords.push_back ({"FILETYPE", {filetype, ""}});
  write_fits (image, file_specific_keywords, file);
}
}
