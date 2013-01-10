#include <vector>
#include <map>

#include <log4cxx/logger.h>
#include <log4cxx/fileappender.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/patternlayout.h>

#include "Params.hxx"
#include "Detector.hxx"
#include "Sample.hxx"
#include "Gnomonic.hxx"
#include "Footprint.hxx"
#include "write_fits.hxx"

log4cxx::LoggerPtr logger(log4cxx::Logger::getRootLogger());

std::map<int,Detector> read_all_DRF_files(const Params::Data_Type &dt,
                                          const std::string &DRF_prefix);

std::vector<Sample> read_all_IN_files(const Params::Data_Type &dt,
                                      const std::string &prefix,
                                      const Gnomonic &projection);

arma::mat make_start_image(const std::string &filename,
                                 const Params &p, int &iter_start);

void compute_correction(const int &nx, const int &ny,
                        const Footprint &fp,
                        const arma::mat &flux_image,
                        const int &iter, const bool &do_cfv,
                        const std::function<double (double)> &boost_func,
                        const int &boost_max_iter,
                        arma::mat &correction,
                        arma::mat &correction_squared);

arma::mat create_spike_image(const int &n, const double &height, const int &nx,
                             const int &ny);

int main(int argc, char* argv[])
{
  Params params(argc,argv);
  Gnomonic projection(params.crval1,params.crval2);

  log4cxx::PatternLayoutPtr layout(new log4cxx::PatternLayout("\%m"));
  log4cxx::FileAppenderPtr
    file_appender(new log4cxx::FileAppender(layout,params.log_filename,false));
  log4cxx::ConsoleAppenderPtr
    console_appender(new log4cxx::ConsoleAppender(layout));
  logger->addAppender(file_appender);
  logger->addAppender(console_appender);

  LOG4CXX_INFO(logger,"HIRES invoked as: ");
  for(int i=0;i<argc;++i)
    LOG4CXX_INFO(logger,argv[i] << " ");
  LOG4CXX_INFO(logger,"\n");

  LOG4CXX_INFO(logger, params);

  std::map<int,Detector> detectors(read_all_DRF_files(params.data_type,
                                                      params.drf_prefix));

  std::vector<Sample> samples(read_all_IN_files(params.data_type,
                                                params.infile_prefix,
                                                projection));
  Footprint footprints(params.radians_per_pix,params.NPIXi,params.NPIXj,
                       params.min_sample_flux,params.angle_tolerance,
                       params.footprints_per_pix,
                       detectors,samples);
  arma::mat wgt_image(footprints.calc_wgt_image(params.NPIXi,
                                                      params.NPIXj));

  if(std::find(params.outfile_types.begin(),params.outfile_types.end(),"cov")
     !=params.outfile_types.end())
    write_fits(wgt_image,"cov",params);

  if(find(params.outfile_types.begin(),params.outfile_types.end(),"flux")
     !=params.outfile_types.end())
    {
      int iter_start;
      arma::mat flux_image=make_start_image(params.starting_image,
                                                  params, iter_start);
      for(int iter= iter_start+1; iter<=params.iter_max; ++iter)
        {
          bool do_cfv_image=
            find(params.outfile_types.begin(),params.outfile_types.end(),"cfv")
            !=params.outfile_types.end()
            && find(params.iter_list.begin(),params.iter_list.end(),iter)
            !=params.iter_list.end();
          arma::mat correction,correction_squared;
          compute_correction(params.NPIXi,params.NPIXj,
                             footprints,flux_image,iter,do_cfv_image,
                             params.boost_func, params.boost_max_iter,
                             correction,correction_squared);
          correction/=wgt_image;
          flux_image%=correction;
          LOG4CXX_INFO(logger,"Mean flux in image " << mean(mean(flux_image))
                       << "\n");
          if(find(params.iter_list.begin(),params.iter_list.end(),iter)
             !=params.iter_list.end())
            write_fits(flux_image,"flux",params,iter);
          if(do_cfv_image)
            {
              arma::mat corr_sq_image = (correction_squared/wgt_image)
                - square(correction);
              write_fits(correction_squared,"cfv",params,iter);
            }
        }
    }
  
  if(find(params.outfile_types.begin(),params.outfile_types.end(),"beam")
     !=params.outfile_types.end())
    {
      arma::mat spike_image=create_spike_image(params.beam_spike_n,
                                               params.beam_spike_height,
                                               params.NPIXi,params.NPIXj);
      // set_fluxes_to_sim_values(all_footprints,spike_image);
      // int iter_start;
      // auto beam_image=make_start_image(params.beam_starting_image,iter_start);
      // for(iter=iter_start+1;iter<=params.iter_max;++iter)
      //   {
      //     Corr_Wgt_Image c(all_footprints,beam_image,iter,false);
      //     beam_image*=c.correction_image;
      //     if(iter_list.find(iter)!=iter_list.end())
      //       write_FITS_image(beam_image,"beam",iter);
      //   }
    }
  LOG4CXX_INFO(logger,"End Processing\n");
}
