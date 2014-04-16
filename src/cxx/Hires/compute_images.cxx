#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
  std::map<int,Detector> read_all_DRF_files(const Hires::Data_Type &dt,
                                            const std::string &DRF_prefix);

  void Hires::compute_images(arma::mat &wgt_image,
                             std::map<int,arma::mat> &flux_images,
                             std::map<int,arma::mat> &cfv_images,
                             std::map<int,arma::mat> &beam_images,
                             std::vector<Sample> &samples,
                             const std::string &outfile_prefix)
  {
    std::map<int,Detector> detectors(read_all_DRF_files(data_type,drf_prefix));

    Footprint footprints(radians_per_pix,ni,nj,min_sample_flux,angle_tolerance,
                         footprints_per_pix,detectors,samples);
    wgt_image=footprints.calc_wgt_image(ni,nj);
    if(std::find(outfile_types.begin(),outfile_types.end(),"cov")
       !=outfile_types.end())
      write_fits(wgt_image,"cov",outfile_prefix);

    if(find(outfile_types.begin(),outfile_types.end(),"flux")
       !=outfile_types.end())
      {
        int iter_start;
        arma::mat flux_image(start_image(starting_image,iter_start));
	arma::mat minimap, hitmap;
     
        if (hires_mode == Hires_Mode::minimap || hires_mode == Hires_Mode::both) {
            footprints.compute_minimap(radians_per_pix,ni,nj,samples,minimap,hitmap);
            write_fits(minimap, "minimap", 0, outfile_prefix);
            write_fits(hitmap, "hitmap", 0, outfile_prefix);
        }

        
        if (hires_mode == Hires_Mode::hires || hires_mode == Hires_Mode::both) {
        for(int iter= iter_start+1; iter<=iter_max; ++iter)
          {
            bool do_cfv_image=
              find(outfile_types.begin(),outfile_types.end(),"cfv")
              !=outfile_types.end();
            arma::mat correction,correction_squared;
            footprints.compute_correction(ni,nj,
                                          flux_image,iter,do_cfv_image,
                                          boost_func, boost_max_iter,
                                          correction,correction_squared);

            correction/=wgt_image;
            flux_image%=correction;

            if(find(iter_list.begin(),iter_list.end(),iter)
               !=iter_list.end())
              {
                flux_images[iter]=flux_image;
                write_fits(flux_image,"flux",iter,outfile_prefix);
                if(do_cfv_image)
                  {
                    arma::mat corr_sq_image = (correction_squared/wgt_image)
                      - square(correction);
                    cfv_images[iter]=corr_sq_image;
                    write_fits(corr_sq_image,"cfv",iter,outfile_prefix);
                  }
              }
          }
        }

      }
  
    if(find(outfile_types.begin(),outfile_types.end(),"beam")
       !=outfile_types.end())
      {
        footprints.set_fluxes_to_sim_values(spike_image());
        int iter_start;
        arma::mat beam_image=(start_image(beam_starting_image,iter_start));
        for(int iter=iter_start+1;iter<=iter_max;++iter)
          {
            arma::mat correction,correction_squared;
            footprints.compute_correction(ni,nj,
                                          beam_image,iter,false,
                                          boost_func, boost_max_iter,
                                          correction,correction_squared);
            beam_image%=correction/wgt_image;
            if(find(iter_list.begin(),iter_list.end(),iter)
               !=iter_list.end())
              {
                beam_images[iter]=beam_image;
                write_fits(beam_image,"beam",iter,outfile_prefix);
              }
          }
      }
  }
}
