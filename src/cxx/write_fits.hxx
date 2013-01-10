#ifndef HIRES_WRITE_FITS_HXX
#define HIRES_WRITE_FITS_HXX

#include <armadillo>
#include <limits>

void write_fits(const arma::mat &image, const std::string &file_type,
                const Params &p, const int &iter);

inline void write_fits(const arma::mat &image, const std::string &file_type,
                       const Params &p)
{
  write_fits(image,file_type,p,std::numeric_limits<int>::max());
}

#endif
