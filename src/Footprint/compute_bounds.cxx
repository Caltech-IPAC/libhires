#include "../Footprint.hxx"

namespace hires
{
  std::vector<int> Footprint::compute_bounds(const arma::mat &response,
                                             const int &i_center,
                                             const int &j_center,
                                             const int &ni, const int &nj)
  {
    int j_size(response.n_rows), i_size(response.n_cols);
    int radius_pixels(j_size/2);
    int j0_resp=0;
    int j1_resp = j_size;
    int i0_resp = 0;
    int i1_resp = i_size;

    /* Compute nominal bounds */
    int j0_image = j_center - radius_pixels;
    int j1_image = j_center + radius_pixels + 1;
    int i0_image = i_center - radius_pixels;
    int i1_image = i_center + radius_pixels + 1;

    /* Trim if outside image */
    if(j0_image<0)
      {
        j0_resp -= j0_image;
        j0_image = 0;
      }
    if(j1_image>nj)
      {
        j1_resp -= j1_image-nj;
        j1_image = nj;
      }
    if(i0_image<0)
      {
        i0_resp -= i0_image;
        i0_image = 0;
      }
    if(i1_image>ni)
      {
        i1_resp -= i1_image-ni;
        i1_image = ni;
      }

    return {j0_image, j1_image, i0_image, i1_image,
        j0_resp, j1_resp, i0_resp, i1_resp};
  }
}
