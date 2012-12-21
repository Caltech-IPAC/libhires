#include "Params.hxx"
// #include "Logfile.hxx"

int main(int argc, char* argv[])
{
  Params params(argc,argv);

  std::cout << params;

  // Logfile log(params.log_filename);

  // log.info << "HIRES invoked as: ";
  // for(int i=0;i<argc;++i)
  //   log.info << argv[i] << " ";
  // log.info << "\n";

  // log.info << params << "\n";
  // const string program(HIRES), version("v1_12");
  // log.step << "Start PROCESSING ("
  //          << program << " "
  //          << version << ")\n";

  // auto all_detectors=read_all_DRF_files(params.DRF_prefix,log);
  // auto all_samples=read_all_IN_files(params.INFILE_prefix,log);
  // auto all_footprints=create_all_footprints(all_samples,all_detectors);
  // auto wgt_image=calc_wgt_image(all_footprints);
  // if(param.outfile_types.find("cov")!=param.outfile_types.end())
  //   write_FITS_image(wgt_image, "cov");

  // if(outfile_types.find("flux")!=outfile_types.end())
  //   {
  //     int iter_start;
  //     auto flux_image=make_start_image(params.starting_image, iter_start);
  //     for(int iter= iter_start+1; iter<=params.iter_max; ++iter)
  //       {
  //         bool do_cfv_image=(params.outfile_types.find("cfv")!=
  //                            params.outfile_types.end())
  //           && params.iter_list.find(iter)!=params.iter_list.end();
  //         Corr_Wgt_Image c(all_footprints,flux_image,iter,do_cfv_image);
  //         flux_image*=c.correction_image;
  //         log.extra <<"Mean flux in image "
  //                   << flux_image.mean()
  //                   << "\n";
  //         if(iter_list.find(iter)!=iter_list.end())
  //           write_FITS_image(flux_image,"flux",iter);
  //         if(do_cfv_image)
  //           {
  //             auto corr_sq_image = (c.sq_wgt / c.wgt) -
  //               (c.correction_image * c.correction_image);
  //             write_FITS_image(corr_sq_image, "cfv", iter);
  //           }
  //       }
  //   }
  
  // if(outfile_types.find("beam")!=outfile_types.end())
  //   {
  //     auto spike_image=create_spike_image(params.beam_spike_n,
  //                                         beam_spike_height);
  //     set_fluxes_to_sim_values(all_footprints,spike_image);
  //     int iter_start;
  //     auto beam_image=make_start_image(params.beam_starting_image,iter_start);
  //     for(iter=iter_start+1;iter<=params.iter_max;++iter)
  //       {
  //         Corr_Wgt_Image c(all_footprints,beam_image,iter,false);
  //         beam_image*=c.correction_image;
  //         if(iter_list.find(iter)!=iter_list.end())
  //           write_FITS_image(beam_image,"beam",iter);
  //       }
  //   }
  // log.step << "End Processing\n";
}
