#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires {

inline std::string GENFILENAME(std::string p, std::string type, int iter) 
{ 
    std::stringstream tmpss;    
    tmpss << p << "_" << type; 
    if (iter != std::numeric_limits<int>::max()) 
        tmpss << "_" << iter;
    tmpss << ".fits"; 
    return (tmpss.str()); 
}

// assumes filename, iter, flux_units from context
#define WRITEFILE(image, desc, isflux) \
{ \
    std::stringstream istr; \
    istr << iter; \
    std::string kwd, val, comment; \
\
    file_specific_keywords.clear(); \
\
    kwd = "ITERNUM"; val = istr.str(); comment = "HIRES iteration number"; \
    file_specific_keywords.push_back( \
        std::make_tuple(kwd, val, comment)); \
\
    kwd = "FILETYPE"; val = desc; comment = ""; \
    file_specific_keywords.push_back( \
        std::make_tuple(kwd, val, comment)); \
\
    if ((isflux)) { \
        kwd = "BUNIT"; val = flux_units; comment = ""; \
        file_specific_keywords.push_back( \
            std::make_tuple(kwd, val, comment)); \
    } \
\
    write_fits((image), file_specific_keywords, filename); \
}


void Hires::write_output(arma::mat &wgt_image,
                             std::map<int,arma::mat> &flux_images,
                             std::map<int,arma::mat> &cfv_images,
                             std::map<int,arma::mat> &beam_images,
                             int &iter,
                             Image_Type image_type,
                             const std::string &outfile_name) 
{
    std::string filename;
    std::vector<std::tuple<std::string,std::string,std::string>> file_specific_keywords;

// Written out only at iteration 0.
    if (iter == 0) {

// Covariance

    if (image_type == Image_Type::cov || 
        (image_type == Image_Type::all &&
            std::find(outfile_types.begin(),outfile_types.end(),"cov")
            != outfile_types.end())) {
        filename = (image_type == Image_Type::all ?
            GENFILENAME(outfile_name, "cov", iter) : 
            outfile_name);
        WRITEFILE(wgt_image, "HIRES covariance image", 0);
    }

// Minimap
    if (image_type == Image_Type::minimap || 
        (image_type == Image_Type::all &&
            std::find(outfile_types.begin(),outfile_types.end(),"minimap")
            != outfile_types.end())) {
        filename = (image_type == Image_Type::all ?
            GENFILENAME(outfile_name, "minimap", iter) : 
            outfile_name);
        WRITEFILE(minimap, "MINIMAP flux image", 1);
    }

// Hit count
    if (image_type == Image_Type::hitmap || 
        (image_type == Image_Type::all &&
            std::find(outfile_types.begin(),outfile_types.end(),"hitmap")
            != outfile_types.end())) {
        filename = (image_type == Image_Type::all ?
            GENFILENAME(outfile_name, "hitmap", iter) : 
            outfile_name);
        WRITEFILE(hitmap, "MINIMAP hit count image", 0);
    }
    } else { 
// Written out at iter > 0
// HIRES flux image
    if (image_type == Image_Type::hires || 
        (image_type == Image_Type::all &&
            std::find(outfile_types.begin(),outfile_types.end(),"hires")
            != outfile_types.end())) {
        filename = (image_type == Image_Type::all ?
            GENFILENAME(outfile_name, "hires", iter) : 
            outfile_name);
        WRITEFILE(flux_images[iter], "HIRES flux image", 1);
    }

// HIRES CFV image
    if (image_type == Image_Type::cfv || 
        (image_type == Image_Type::all &&
            std::find(outfile_types.begin(),outfile_types.end(),"cfv")
            != outfile_types.end())) {
        filename = (image_type == Image_Type::all ?
            GENFILENAME(outfile_name, "cfv", iter) : 
            outfile_name);
        WRITEFILE(cfv_images[iter], "HIRES correction factor variance image", 0);
    }

// HIRES beam image
    if (image_type == Image_Type::beam || 
        (image_type == Image_Type::all &&
            std::find(outfile_types.begin(),outfile_types.end(),"beam")
            != outfile_types.end())) {
        filename = (image_type == Image_Type::all ?
            GENFILENAME(outfile_name, "beam", iter) : 
            outfile_name);
        WRITEFILE(beam_images[iter], "HIRES beam image", 0);
    }
    } // end iter > 0
}
} // end namespace hires

