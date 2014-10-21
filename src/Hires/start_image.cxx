#include <CCfits/CCfits>
#include <valarray>
#include "../Hires.hxx"
#include "../Exception.hxx"

namespace hires
{
Eigen::MatrixXd Hires::start_image (const std::string &filename, int &iter_start)
{
  Eigen::MatrixXd image (nxy[1], nxy[0]);
  iter_start = 0;
  if (filename.empty())
    {
      image.fill (1.0);
    }
  else
    {
      CCfits::FITS hdus (filename);
      CCfits::PHDU &phdu (hdus.pHDU ());
      phdu.readAllKeys ();
      if (nxy[0] != phdu.axis (0) || nxy[1] != phdu.axis (1))
        {
          std::stringstream ss;
          ss << "STARTING_IMAGE " << filename << "has incorrect dimensions\n"
             << "  Must be: " << nxy[0] << " " << nxy[1] << "\n"
             << "  Actual value: " << phdu.axis (0) << " " << phdu.axis (1);
          throw Exception (ss.str ());
        }
      std::map<std::string, CCfits::Keyword *> &m (phdu.keyWord ());
      double CRVAL1, CRVAL2;
      m["CRVAL1"]->value (CRVAL1);
      if (std::abs (CRVAL1 - crval[0]) > 0.001)
        {
          std::stringstream ss;
          ss << "STARTING_IMAGE " << filename << " has inconsistent values\n"
             << "  Should be: " << crval[0] << "\n"
             << "  Actual value: " << CRVAL1;
          throw Exception (ss.str ());
        }
      m["CRVAL2"]->value (CRVAL2);
      if (std::abs (CRVAL2 - crval[1]) > 0.001)
        {
          std::stringstream ss;
          ss << "STARTING_IMAGE " << filename << " has inconsistent CRVAL1\n"
             << "  Should be: " << crval[1] << "\n"
             << "  Actual value: " << CRVAL2;
          throw Exception (ss.str ());
        }

      if (m.find ("ITERNUM") != m.end ())
        m["ITERNUM"]->value (iter_start);

      std::valarray<double> valarray_image;
      phdu.read (valarray_image);
      for (int i = 0; i < nxy[0]; ++i)
        for (int j = 0; j < nxy[1]; ++j)
          image (j, i) = valarray_image[i + nxy[0] * j];
    }
  return image;
}
}
