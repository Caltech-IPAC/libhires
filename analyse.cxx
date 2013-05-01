#include <CCfits/CCfits>

int main(int argc, char *argv[])
{
  if(argc!=3)
    {
      std::cerr << "Need exactly two arguments\n";
      exit(1);
    }
  CCfits::FITS py(argv[1]);
  CCfits::FITS cxx(argv[2]);

  std::valarray<double> py_image,cxx_image;
  py.pHDU().read(py_image);
  cxx.pHDU().read(cxx_image);

  std::valarray<double> diff(py_image-cxx_image);

  std::cout << "Sum: " << std::abs(py_image+cxx_image).max()/2 << "\n"
            << "Diff: " << std::abs(diff).max() << "\n"
            << "Scaled Diff: "
            << std::abs(diff/(py_image+cxx_image+abs(cxx_image).max()*1e-30)).max() << "\n";
}
