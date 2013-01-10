#include <armadillo>
#include <CCfits>
#include <valarray>
#include "Params.hxx"
#include "logger.hxx"

arma::mat make_start_image(const std::string &filename,
                                 const Params &p, int &iter_start)
{
  arma::mat image(p.NPIXi,p.NPIXj);
  iter_start=0;
  if(filename=="flat")
    {
      image.fill(1.0);
    }
  else
    {
      CCfits::FITS hdus(filename);
      CCfits::PHDU &phdu(hdus.pHDU());
      phdu.readAllKeys();
      if(p.NPIXi!=phdu.axis(0) || p.NPIXj!=phdu.axis(1))
        LOG4CXX_FATAL(logger,"STARTING_IMAGE " << filename
                      << "has incorrect dimensions\n"
                      << "  Must be: " << p.NPIXi << " " << p.NPIXj << "\n"
                      << "  Actual value: " << phdu.axis(0) << " "
                      << phdu.axis(1) << "\n");
      std::map<std::string,CCfits::Keyword *> &m(phdu.keyWord());
      double crval1, crval2;
      std::string bunit;
      m["CRVAL1"]->value(crval1);
      if(std::abs(crval1-p.crval1)>0.001)
        LOG4CXX_FATAL(logger,"STARTING_IMAGE " << filename
                      << " has inconsistent values\n"
                      << "  Should be: " << p.crval1 << "\n"
                      << "  Actual value: " << crval1 << "\n");
      m["CRVAL2"]->value(crval2);
      if(std::abs(crval2-p.crval2)>0.001)
        LOG4CXX_FATAL(logger,"STARTING_IMAGE " << filename
                      << " has inconsistent CRVAL1\n"
                      << "  Should be: " << p.crval2 << "\n"
                      << "  Actual value: " << crval2 << "\n");
      m["BUNIT"]->value(bunit);
      if(bunit!=p.flux_units)
        LOG4CXX_FATAL(logger,"STARTING_IMAGE " << filename
                      << " has inconsistent BUNIT\n"
                      << "  Should be: " << p.flux_units << "\n"
                      << "  Actual value: " << bunit << "\n");
      if(m.find("ITERNUM")!=m.end())
        {
          m["ITERNUM"]->value(iter_start);
          LOG4CXX_INFO(logger,"ITERNUM in " << filename << " is "
                       << iter_start << "\n");
        }
      else
        {
          LOG4CXX_INFO(logger,"No ITERNUM keyword in START_IMAGE "
                       << filename << "\n");
        }
        
      std::valarray<double> valarray_image;
      phdu.read(valarray_image);
      for(int i=0;i<p.NPIXi;++i)
        for(int j=0;j<p.NPIXj;++j)
          image(i,j)=valarray_image[i+p.NPIXi*j];
    }
  return image;
}
