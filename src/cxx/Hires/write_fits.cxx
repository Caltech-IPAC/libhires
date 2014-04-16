/* Writes out one FITS image of specified type.
   image = array with pixel values
   file_type = 'flux', 'cov', 'cfv', etc.
   iter = iteration number (put in FITS keyword)
*/

#include <chrono>
#include <limits>
#include <CCfits/CCfits>
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include "../Hires.hxx"
#include "../version.hxx"

namespace hires
{
  void Hires::write_fits(const arma::mat &image, const std::string &file_type,
                         const int &iteration, const std::string &outfile_prefix)
  {
    if(outfile_prefix.empty())
      return;
    boost::filesystem::path dir(outfile_prefix);
    boost::filesystem::create_directory(dir);
    std::stringstream filename;
    filename << outfile_prefix << "_" << file_type;
    if(iteration!=std::numeric_limits<int>::max())
      filename << "_" << iteration;
    filename << ".fits";
    boost::filesystem::path fits_file(filename.str());
  
    long axes[]={image.n_cols,image.n_rows};
    boost::filesystem::remove(fits_file);
    CCfits::FITS outfile(fits_file.string(),FLOAT_IMG,2,axes);

    CCfits::PHDU &phdu(outfile.pHDU());

    if(file_type=="flux")
      phdu.addKey("BUNIT",flux_units,"");
    for(auto &keywords: fits_keywords)
      phdu.addKey(std::get<0>(keywords),std::get<1>(keywords),
                  std::get<2>(keywords));
    phdu.addKey("CRVAL1",crval1,"");
    phdu.addKey("CRVAL2",crval2,"");
    phdu.addKey("CTYPE1",ctype1,"");
    phdu.addKey("CTYPE2",ctype2,"");
    float cdelt_rounded(radians_per_pix*180/boost::math::constants::pi<double>());
    phdu.addKey("CD1_1",-cdelt_rounded,"Degrees / Pixel");
    phdu.addKey("CD2_1",0.0,"Degrees / Pixel");
    phdu.addKey("CD1_2", -0.0, "Degrees / Pixel");
    phdu.addKey("CD2_2",cdelt_rounded,"Degrees / Pixel");

    phdu.addKey("CRPIX1", (ni+1)/2 , "center pixel");
    phdu.addKey("CRPIX2", (nj+1)/2 , "center pixel");
    std::string t_comment("HIRES");
    if(file_type=="flux")
      t_comment+=" flux";
    else if(file_type=="cov")
      t_comment+= " coverage";
    else if(file_type=="cfv")
      t_comment+= " correction factor variance";
    t_comment+= " image";
    phdu.addKey("FILETYPE", file_type, t_comment);
    if(iteration!=std::numeric_limits<int>::max())
      phdu.addKey("ITERNUM", iteration, "HIRES iteration number");
    phdu.addKey("FILENAME",
                fits_file.filename().string(),
                "name of this file");
    if(!drf_prefix.empty())
      phdu.addKey("DRF_IN",
                  boost::filesystem::path(drf_prefix).filename().string()+"*",
                  "Name of Detector Response Files");

    std::chrono::system_clock::time_point tp(std::chrono::system_clock::now());
    std::time_t now(std::chrono::system_clock::to_time_t(tp));
    phdu.addKey("DATE",std::string(std::ctime(&now)),
                "when this file was created");
    phdu.addKey("CREATED","HIRES " + version,
                "software version that created this file");

    std::valarray<float> temp(image.n_cols*image.n_rows);
    for(size_t i=0;i<image.n_cols;++i)
      for(size_t j=0;j<image.n_rows;++j)
// GLON runs backwards
        temp[i+image.n_cols*j]=image(j,(image.n_cols - 1) - i);
    phdu.write(1,temp.size(),temp);
  }
}
