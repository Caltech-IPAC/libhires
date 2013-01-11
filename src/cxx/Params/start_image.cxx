#include <armadillo>
#include <CCfits>
#include <valarray>
#include "../Params.hxx"
#include "../Exception.hxx"

namespace hires
{
  arma::mat Params::start_image(const std::string &filename, int &iter_start)
  {
    arma::mat image(nj,ni);
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
        if(ni!=phdu.axis(0) || nj!=phdu.axis(1))
          {
            std::stringstream ss;
            ss << "STARTING_IMAGE " << filename
               << "has incorrect dimensions\n"
               << "  Must be: " << ni << " " << nj << "\n"
               << "  Actual value: " << phdu.axis(0) << " "
               << phdu.axis(1);
            throw Exception(ss.str());
          }
        std::map<std::string,CCfits::Keyword *> &m(phdu.keyWord());
        double CRVAL1, CRVAL2;
        m["CRVAL1"]->value(CRVAL1);
        if(std::abs(CRVAL1-crval1)>0.001)
          {
            std::stringstream ss;
            ss << "STARTING_IMAGE " << filename
               << " has inconsistent values\n"
               << "  Should be: " << crval1 << "\n"
               << "  Actual value: " << CRVAL1;
            throw Exception(ss.str());
          }
        m["CRVAL2"]->value(CRVAL2);
        if(std::abs(CRVAL2-crval2)>0.001)
          {
            std::stringstream ss;
            ss << "STARTING_IMAGE " << filename
               << " has inconsistent CRVAL1\n"
               << "  Should be: " << crval2 << "\n"
               << "  Actual value: " << CRVAL2;
            throw Exception(ss.str());
          }

        if(m.find("ITERNUM")!=m.end())
          m["ITERNUM"]->value(iter_start);
        
        std::valarray<double> valarray_image;
        phdu.read(valarray_image);
        for(int i=0;i<ni;++i)
          for(int j=0;j<nj;++j)
            image(j,i)=valarray_image[i+ni*j];
      }
    return image;
  }
}
