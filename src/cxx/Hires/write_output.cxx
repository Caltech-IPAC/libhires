#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires {
void Hires::write_output(arma::mat &wgt_image,
                             std::map<int,arma::mat> &flux_images,
                             std::map<int,arma::mat> &cfv_images,
                             std::map<int,arma::mat> &beam_images,
                             int &iter,
                             const std::string &outfile_prefix) 
{
    if (iter == 0 &&
        std::find(outfile_types.begin(),outfile_types.end(),"cov")
       != outfile_types.end())
        write_fits(wgt_image,"cov",outfile_prefix);

    if (find(outfile_types.begin(),outfile_types.end(),"flux")
       != outfile_types.end()) {
     
        if (iter == 0 &&
            (hires_mode == Hires_Mode::minimap || hires_mode == Hires_Mode::both)) {
            write_fits(minimap, "minimap", 0, outfile_prefix);
            write_fits(hitmap, "hitmap", 0, outfile_prefix);
        }
        
        if (iter != 0 && 
            (hires_mode == Hires_Mode::hires || hires_mode == Hires_Mode::both)) {
            write_fits(flux_images[iter],"flux",iter,outfile_prefix);
            if (find(outfile_types.begin(),outfile_types.end(),"cfv") 
                != outfile_types.end())
                write_fits(cfv_images[iter],"cfv",iter,outfile_prefix);
        }
    }
  
    if (iter != 0 &&
        find(outfile_types.begin(),outfile_types.end(),"beam")
            != outfile_types.end())
        write_fits(beam_images[iter],"beam",iter,outfile_prefix);
}
} // end namespace Hires
