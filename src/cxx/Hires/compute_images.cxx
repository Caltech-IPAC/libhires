#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
void Hires::compute_images(arma::mat &wgt_image,
                             std::map<int,arma::mat> &flux_images,
                             std::map<int,arma::mat> &cfv_images,
                             std::map<int,arma::mat> &beam_images,
                             std::vector<Sample> &samples,
                             const std::string &outfile_prefix) 
{
    for (int iter=0; iter<=iter_max; ++iter) {
        iterate(wgt_image, flux_images, cfv_images, beam_images, 
                samples, iter);

        if (iter == 0 || find(iter_list.begin(),iter_list.end(),iter) != iter_list.end()) 
            write_output(wgt_image, flux_images, cfv_images, beam_images, 
                         iter, outfile_prefix);
    }
}
} // end namespace hires
