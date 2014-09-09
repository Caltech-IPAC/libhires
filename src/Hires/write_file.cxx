#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{

void Hires::write_file (arma::mat image, std::string filename,
                        const char *desc, int iter, int isflux)
{
  std::stringstream istr;
  istr << iter;
  std::string kwd, val, comment;
  std::vector<std::tuple<std::string, std::string, std::string> >
  file_specific_keywords;

  file_specific_keywords.clear ();

  kwd = "ITERNUM";
  val = istr.str ();
  comment = "HIRES iteration number";
  file_specific_keywords.push_back (std::make_tuple (kwd, val, comment));

  kwd = "FILETYPE";
  val = desc;
  comment = "";
  file_specific_keywords.push_back (std::make_tuple (kwd, val, comment));

  if (iter!=0)
    file_specific_keywords.push_back
      (std::make_tuple ("DRF_IN",
                        boost::filesystem::path (drf_prefix).filename ().string ()
                        + "*",
                        "Name of Detector Response Files"));
  if (isflux)
    {
      kwd = "BUNIT";
      val = flux_units;
      comment = "";
      file_specific_keywords.push_back (std::make_tuple (kwd, val, comment));
    }

  write_fits (image, file_specific_keywords, filename);
}
} // end namespace hires
