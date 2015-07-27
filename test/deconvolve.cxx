#include <random>
#include <boost/filesystem.hpp>
#include <iostream>
#include <CCfits/CCfits>

#include "../src/Hires.hxx"

bool check_exists(const boost::filesystem::path &output)
{
  if (!exists(output))
    {
      std::cout << "FAIL: " << output.string() << " does not exist\n";
      return false;
    }
  return true;
}
  
bool compare_fits(const boost::filesystem::path &computed_name,
                  const boost::filesystem::path &expected_name)
{
  CCfits::FITS computed(computed_name.c_str());
  CCfits::FITS expected(expected_name.c_str());

  std::valarray<double> computed_image,expected_image;
  computed.pHDU().read(computed_image);
  expected.pHDU().read(expected_image);

  std::valarray<double> diff(computed_image-expected_image);
  bool passed=(std::abs(diff/(computed_image + expected_image
                              + abs(expected_image).max()*1e-30)).max() < 1e-14);
  if (passed)
    std::cout << "PASS: " << expected_name.stem() << "\n";
  else
    std::cout << "FAIL: " << expected_name.stem() << "\n";
  return passed;
}

int main()
{
  /// Generate a bunch of pseudo-random samples and run hires.
  std::mt19937 gen(0);
  const std::array<size_t,2> nxy{{32,32}};
  const std::array<double,2> crval{{0,0}};
  const double radians_per_pix(0.001);
  std::uniform_real_distribution<> uniform(-(nxy[0]*radians_per_pix/2),
                                           nxy[0]*radians_per_pix/2);
  const double x0(7*radians_per_pix), y0(10*radians_per_pix),
    sigma_x2(3.14*3.14*radians_per_pix*radians_per_pix),
    sigma_y2(2.72*2.72*radians_per_pix*radians_per_pix),
    noise_rms(0.01);
  std::normal_distribution<> noise(0,noise_rms);
  std::vector<hires::Sample> samples;
  const size_t num_samples(10000);
  for (size_t s=0; s<num_samples; ++s)
    {
      double x(uniform(gen)), y(uniform(gen));
      samples.emplace_back(x,y,
                           std::exp(-((x-x0)*(x-x0)/sigma_x2
                                      + (y-y0)*(y-y0)/sigma_y2))
                           + noise(gen),0);
    }
  hires::Hires hires (nxy,crval,0.001,std::vector<std::pair<std::string,
                      std::pair<std::string,std::string> > >(),samples);
  hires.compute_minimap();
  const double sigma_drf(3.5*radians_per_pix);
  const int num_iterations(2);
  hires.compute_mcm(sigma_drf,num_iterations);
  hires.compute_elastic_net(sigma_drf);
  
  boost::filesystem::path minimap("test/minimap/minimap.fits"),
    mcm("test/mcm/mcm.fits"),
    elastic_net("test/elastic_net/elastic_net.fits");
  hires.write_minimap("test/minimap/");
  hires.write_mcm("test/mcm/");
  hires.write_elastic_net("test/elastic_net/");
  if (!check_exists(minimap)
      || !check_exists(mcm)
      || !check_exists(elastic_net))
    exit(-1);

  bool result=compare_fits(minimap,"test/expected/minimap.fits");
  result=compare_fits(mcm,"test/expected/mcm.fits") && result;
  result=compare_fits(elastic_net,"test/expected/elastic_net.fits") && result;

  remove(minimap);
  remove(minimap.parent_path());
  remove(mcm);
  remove(mcm.parent_path());
  remove(elastic_net);
  remove(elastic_net.parent_path());
  return result ? 0 : 1;
}
