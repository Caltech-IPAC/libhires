#include <thread>

#include "../Footprint.hxx"

namespace hires
{
void Footprint::compute_correction (const Eigen::MatrixXd &signal_image,
                                    Eigen::MatrixXd &correction)
{
  const Eigen::MatrixXd &AA=response;

  Eigen::MatrixXd ddB(result.size(),result.size());
  ddB.setZero();
  Eigen::VectorXd dB(result.size());

  double result_norm=result.norm(), du_norm=0;

  Eigen::VectorXd dA(2*AA*result - rhs);

  Eigen::VectorXd entropy(result.size());

  for(int i=0; i<result.size(); ++i)
    {
      entropy(i)=log(std::cosh(result(i)/noise_level));
      dB(i)=std::tanh(result(i)/noise_level)/noise_level;
      double temp=std::cosh(result(i)/noise_level)*noise_level;
      ddB(i,i)=1/(temp*temp);
    }

  Eigen::VectorXd du=-(AA + lambda*ddB).fullPivHouseholderQr()
    .solve(dA + lambda*dB);

  // FIXME: Doing the correction in two places.
  result+=du;

  for(int i=0; i<correction.rows(); ++i)
    for(int j=0; j<correction.cols(); ++j)
      correction(i,j)=du(i+correction.rows()*j);

  result_norm=result.norm();
  du_norm=du.norm();

  std::cout << result_norm << " "
            << du_norm << "\n";
}
}
