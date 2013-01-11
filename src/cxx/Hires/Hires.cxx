#include "../Hires.hxx"
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include "../Exception.hxx"

namespace hires
{
  Hires::Hires(const std::string &Data_type, const std::string &Infile_prefix,
               const std::string &Outfile_prefix,
               const std::vector<std::string> &param_files):
    infile_prefix(Infile_prefix),
    outfile_prefix(Outfile_prefix),
    starting_image("flat"),
    beam_starting_image("flat"),
    flux_units("??"),
    log_filename("logfile.log"),
    outfile_types{"flux"},
    boost_max_iter(0),
    footprints_per_pix(1),
    beam_spike_n(5),
    min_sample_flux(std::numeric_limits<double>::min()),
    beam_spike_height(10)
    {
      std::map<std::string,Data_Type> data_types{{"planck",Data_Type::planck},
          {"spire",Data_Type::spire}};
      if(data_types.find(Data_type)!=data_types.end())
        data_type=data_types[Data_type];
      else
        {
          std::stringstream ss;
          ss << "data_type is '" << Data_type
             << "' but must be one of ( ";
          for(auto &m: data_types)
            ss << "'" << m.first << "' ";
          ss << ")";
          throw Exception(ss.str());
        }

      if(data_type==Data_Type::planck)
        flux_units="K_RJ";
      else if(data_type==Data_Type::spire)
        flux_units="Jy/beam";

      for(auto &pf: param_files)
        {
          std::ifstream infile(pf);
          while(infile)
            {
              std::string s;
              getline(infile,s);
              std::vector<std::string> line;
              boost::split(line,s,boost::is_any_of("#"));
              if(line.empty() || line[0].empty())
                continue;
              std::vector<std::string> words;
              boost::split(words,line[0],boost::is_any_of("\t "),
                           boost::token_compress_on);
              if(words.size()<2)
                throw Exception("Parameter: " + words[0]
                                + " has no value in file: " + pf + "\n");
              std::stringstream ss(line[0]);
              std::string key;
              ss >> key;
              if(key=="SIZE_NPIX")
                {
                  ss >> ni;
                  if(words.size()>=3)
                    ss >> nj;
                  else
                    nj=ni;
                }
              else if(key=="ARCSEC_PER_PIX")
                {
                  ss >> radians_per_pix;
                  radians_per_pix*=boost::math::constants::pi<double>()
                    /(3600*180);
                }
              else if(key=="CRVAL1")
                ss >> crval1;
              else if(key=="CRVAL2")
                ss >> crval2;
              else if(key=="CTYPE1")
                ss >> ctype1;
              else if(key=="CTYPE2")
                ss >> ctype2;
              else if(key=="OUTFILE_TYPES")
                outfile_types=std::vector<std::string>(words.begin()+1,
                                                       words.end());
              else if(key=="ANGLE_TOLERANCE")
                {
                  ss >> angle_tolerance;
                  angle_tolerance*=boost::math::constants::pi<double>()/180;
                }
              else if(key=="FOOTPRINTS_PER_PIX")
                ss >> footprints_per_pix;
              else if(key=="DRF_PREFIX")
                ss >> drf_prefix;
              else if(key=="STARTING_IMAGE")
                ss >> starting_image;
              else if(key=="MIN_SAMPLE_FLUX")
                ss >> min_sample_flux;
              else if(key=="BEAM_SPIKE_N")
                ss >> beam_spike_n;
              else if(key=="BEAM_SPIKE_HEIGHT")
                ss >> beam_spike_height;
              else if(key=="BEAM_STARTING_IMAGE")
                ss >> beam_starting_image;
              else if(key=="FLUX_UNITS")
                ss >> flux_units;
              else if(key=="LOG_FILENAME")
                ss >> log_filename;
              else if(key=="ITER_MAX")
                ss >> iter_max;
              else if(key=="ITER_LIST")
                {
                  for(auto w=words.begin()+1; w!=words.end(); ++w)
                    {
                      std::vector<std::string> iter_str;
                      boost::split(iter_str,*w,boost::is_any_of(","));
                      for(auto &n_str: iter_str)
                        {
                          if(!n_str.empty())
                            {
                              std::stringstream ns(n_str);
                              int n;
                              ns >> n;
                              iter_list.push_back(n);
                            }
                        }
                    }
                }
              else if(key=="BOOST_CORRECTION")
                {
                  ss >> boost_max_iter;
                  if(words.size() >= 3)
                    ss >> boost_type;
                  if(boost_type == "TIMES_2")
                    boost_func = [](const double &x) {return x+x-1.0;};
                  else if(boost_type == "TIMES_3")
                    boost_func = [](const double &x) {return x+x+x-2.0;};
                  else if(boost_type == "SQUARED")
                    boost_func = [](const double &x) {return x*x;};
                  else if(boost_type == "EXP_2.5")
                    boost_func = [](const double &x) {return pow(x, 2.5);};
                  else if(boost_type == "CUBED")
                    boost_func = [](const double &x) {return x*x*x;};
                  else
                    throw Exception("Unknown BOOST type: " + boost_type);
                }
              else if(key=="KWD")
                {
                  std::string kwd;
                  ss >> kwd;
                  bool need_warn(true);

                  if(kwd=="CRVAL1")
                    ss >> crval1;
                  else if(kwd == "CRVAL2")
                    ss >> crval2;
                  else if(kwd == "CTYPE1")
                    ss >> ctype1;
                  else if(kwd == "CTYPE2")
                    ss >> ctype2;
                  else
                    {
                      std::string val;
                      ss >> val;
                      std::string comment;
                      getline(ss,comment);
                      fits_keywords.push_back(std::make_tuple(kwd,val,comment));
                      need_warn = false;
                    }
                  if(need_warn)
                    {
                      std::cerr << "Using KWD before " << kwd
                                << " parameter is deprecated\n";
                    }
                }
              else if(key=="")
                continue;
              else
                throw Exception("ERROR: Unknown parameter name: "
                                + key + " in file: " + pf);
            }
        }
      /* check for errors */
      std::vector<std::string> valid_outfile_types={"flux","cov","beam","cfv"};
      for(auto &t: outfile_types)
        {
          if(find(valid_outfile_types.begin(),valid_outfile_types.end(),t)
             ==valid_outfile_types.end())
            throw Exception("Illegal OUTFILE_TYPE: " + t);
        }

      iter_list.push_back(iter_max); /* make sure to output files for
                                        final iteration */
      std::sort(iter_list.begin(),iter_list.end());
      iter_list.resize(std::unique(iter_list.begin(),
                                   iter_list.end())-iter_list.begin());
    }
}