#include "../Hires.hxx"
#include "../Detector.hxx"
#include "../Footprint.hxx"

namespace hires
{
std::map<int,Detector> read_all_DRF_files(const std::string &DRF_prefix);

void Hires::iterate(arma::mat &wgt_image,
                             std::map<int,arma::mat> &flux_images,
                             std::map<int,arma::mat> &cfv_images,
                             std::map<int,arma::mat> &beam_images,
                             std::vector<Sample> &samples,
			     int &iter) 
{
    std::map<int,Detector> detectors(read_all_DRF_files(drf_prefix));

    Footprint footprints(radians_per_pix,ni,nj,min_sample_flux,angle_tolerance,
                         footprints_per_pix,detectors,samples);
    wgt_image = footprints.calc_wgt_image(ni,nj);
    int iter_start;

//
// Standard hires/minimap processing - always happens.
//
    if (iter == 0) {
        flux_images[0] = start_image(starting_image,iter_start);
        footprints.compute_minimap(radians_per_pix,ni,nj,samples,minimap,hitmap);
    } else {
        arma::mat correction,correction_squared;
        footprints.compute_correction(ni,nj,
                                      flux_images[iter-1],iter, 1,
                                      boost_func, boost_max_iter,
                                      correction,correction_squared);
        correction /= wgt_image;
        flux_images[iter] = flux_images[iter-1] % correction; // Schur product

        arma::mat corr_sq_image = (correction_squared/wgt_image) - square(correction);
        cfv_images[iter]=corr_sq_image;
    }

//
// Optional HIRES beam generation
//
    if (find(outfile_types.begin(),outfile_types.end(),"beam")
        != outfile_types.end()) {
          footprints.set_fluxes_to_sim_values(spike_image());

        if (iter == 0) {
            beam_images[0] = start_image(beam_starting_image,iter_start);
        } else {
            arma::mat correction,correction_squared;
            footprints.compute_correction(ni,nj,
                                      beam_images[iter-1],iter,false,
                                      boost_func, boost_max_iter,
                                      correction,correction_squared);
            beam_images[iter] = beam_images[iter-1] % correction/wgt_image; // Schur product
        }
    }
}
} // end namespace hires
