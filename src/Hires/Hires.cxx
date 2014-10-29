#include "../Hires.hxx"
#include "../Footprint.hxx"

namespace hires
{
Hires::Hires (const std::array<int,2> &Nxy,
              const std::array<double,2> &Crval, const double &Radians_per_pix,
              const std::set<Image_Type> &Output_types,
              const boost::filesystem::path &Drf_file,
              const std::vector<std::pair<std::string, std::pair<std::string,
                                                                 std::string> > >
              &Fits_keywords,
              const std::vector<Sample> &Samples):
  drf_file(Drf_file),
  nxy(Nxy), footprints_per_pix (1),
  crval(Crval), radians_per_pix(Radians_per_pix),
  angle_tolerance (2.5),
  output_types(Output_types),
  fits_keywords(Fits_keywords),
  samples(Samples),
  iteration(0)
{
  fits_keywords.push_back({std::string ("CREATOR"),
        {std::string ("LIBHIRES"), std::string("")}});
              
  if(output_types.find(Image_Type::minimap_image)!=output_types.end()
     || output_types.find(Image_Type::minimap_hitmap)!=output_types.end())
    compute_minimap ();

  if(running_hires())
    {
      detectors=read_DRF (drf_file);
      footprints=Footprint (radians_per_pix, nxy, angle_tolerance,
                            footprints_per_pix, detectors, samples);
      signal_image.resize(nxy[0],nxy[1]);
      signal_image.setZero();
    }
}
}
