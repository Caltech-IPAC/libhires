#include "../Footprint.hxx"

namespace hires
{
std::vector<int> Footprint::compute_bounds (const Eigen::MatrixXd &response,
                                            const int &i_center,
                                            const int &j_center,
                                            const std::array<int,2> &nxy) const
{
  // FIXME: This assumes that num_columns==num_rows
  int j_size (response.rows()), i_size (response.cols());
  int radius_pixels (j_size / 2);
  int j0_resp = 0;
  int j1_resp = j_size;
  int i0_resp = 0;
  int i1_resp = i_size;

  /* Compute nominal bounds */
  int j0_image = j_center - radius_pixels;
  // FIXME: Is this correct?  It only works if response.rows() is odd.
  int j1_image = j_center + radius_pixels + 1;
  int i0_image = i_center - radius_pixels;
  int i1_image = i_center + radius_pixels + 1;

  /* Trim if outside image */
  if (j0_image < 0)
    {
      j0_resp -= j0_image;
      j0_image = 0;
    }
  if (j1_image > nxy[1])
    {
      j1_resp -= j1_image - nxy[1];
      j1_image = nxy[1];
    }
  if (i0_image < 0)
    {
      i0_resp -= i0_image;
      i0_image = 0;
    }
  if (i1_image > nxy[0])
    {
      i1_resp -= i1_image - nxy[0];
      i1_image = nxy[0];
    }

  return { j0_image, j1_image, i0_image, i1_image,
           j0_resp,  j1_resp,  i0_resp,  i1_resp };
}
}
