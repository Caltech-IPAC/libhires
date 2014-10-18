#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{

inline std::string genfilename (std::string p, std::string type, int iter)
{
  std::stringstream tmpss;
  tmpss << p << "_" << type;
  if (iter != std::numeric_limits<int>::max ())
    tmpss << "_" << iter;
  tmpss << ".fits";
  return (tmpss.str ());
}

void Hires::write_output (Image_Type image_type,
                          const std::string &outfile_name)
{
  std::string filename;

  if (iteration == 0)
    {
      // Covariance
      if (image_type == Image_Type::cov || image_type == Image_Type::all)
        {
          filename = (image_type == Image_Type::all
                          ? genfilename (outfile_name, "cov", iteration)
                          : outfile_name);
          write_file (wgt_image, filename, "HIRES covariance image", iteration, 0);
        }

      // Minimap
      if (image_type == Image_Type::minimap || image_type == Image_Type::all)
        {
          filename = (image_type == Image_Type::all
                          ? genfilename (outfile_name, "minimap", iteration)
                          : outfile_name);
          write_file (minimap, filename, "MINIMAP flux image", iteration, 1);
        }

      // Hit count
      if (image_type == Image_Type::hitmap || image_type == Image_Type::all)
        {
          filename = (image_type == Image_Type::all
                          ? genfilename (outfile_name, "hitmap", iteration)
                          : outfile_name);
          write_file (hitmap, filename, "MINIMAP hit count image", iteration, 0);
        }
    }
  else
    {
      // Written out at iter > 0
      // HIRES flux image
      if (image_type == Image_Type::hires || image_type == Image_Type::all)
        {
          filename = (image_type == Image_Type::all
                          ? genfilename (outfile_name, "hires", iteration)
                          : outfile_name);
          write_file (flux_images, filename, "HIRES flux image", iteration,
                      1);
        }

      // HIRES CFV image
      if (image_type == Image_Type::cfv || image_type == Image_Type::all)
        {
          filename = (image_type == Image_Type::all
                          ? genfilename (outfile_name, "cfv", iteration)
                          : outfile_name);
          write_file (cfv_images, filename,
                      "HIRES correction factor variance image", iteration, 0);
        }

      // HIRES beam image
      if (image_type == Image_Type::beam
          || (image_type == Image_Type::all && generate_beams))
        {
          if(!generate_beams)
            throw Exception("Asking for beam output, but no beams were "
                            "generated");
          filename = (image_type == Image_Type::all
                          ? genfilename (outfile_name, "beam", iteration)
                          : outfile_name);
          write_file (beam_images, filename, "HIRES beam image", iteration,
                      0);
        }
    }
}
}
