#include "../Hires.hxx"
#include "../Footprint.hxx"

namespace hires
{
Hires::Hires (const std::array<int,2> &Nxy,
              const std::array<double,2> &Crval, const double &Radians_per_pix,
              const std::set<Image_Type> &Output_types,
              const boost::filesystem::path &Drf_file,
              const std::string &boost_function_string,
              const std::vector<std::pair<std::string, std::pair<std::string,
                                                                 std::string> > >
              &Fits_keywords,
              const std::vector<Sample> &Samples):
  drf_file(Drf_file),
  nxy(Nxy), footprints_per_pix (1), beam_spike_n (5),
  crval(Crval), radians_per_pix(Radians_per_pix),
  angle_tolerance (2.5),
  beam_spike_height (10),
  output_types(Output_types),
  fits_keywords(Fits_keywords),
  samples(Samples),
  iteration(0)
{
  fits_keywords.push_back({std::string ("CREATOR"),
        {std::string ("LIBHIRES"), std::string("")}});
              

  std::map<std::string,std::function<double(double)> > boost_functions=
    {{"TIMES_2",[](const double &x) { return x + x - 1.0; }},
     {"TIMES_3", [](const double &x) { return x + x + x - 2.0; }},
     {"SQUARED", [](const double &x) { return x * x; }},
     {"EXP_2.5", [](const double &x) { return pow (x, 2.5); }},
     {"CUBED", [](const double &x) { return x * x * x; }}};

  auto f=boost_functions.find(boost_function_string);
  if(f==boost_functions.end())
    {
      if(!boost_function_string.empty())
        throw Exception("Invalid boost function string: "
                        + boost_function_string);
    }
  else
    {
      boost_function=f->second;
    }

  if(output_types.find(Image_Type::minimap_image)!=output_types.end()
     || output_types.find(Image_Type::minimap_hitmap)!=output_types.end())
    compute_minimap ();

  if(running_hires())
    {
      detectors=read_DRF (drf_file);
      footprints=Footprint (radians_per_pix, nxy, angle_tolerance,
                            footprints_per_pix, detectors, samples);
      weight_image = footprints.calc_wgt_image (nxy);
      // FIXME: Shouldn't I just use iteration?
      int iter_start;

      signal_image = start_image (starting_image, iter_start);
    }
}
}
