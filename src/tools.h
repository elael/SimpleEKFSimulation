#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {
 public:
  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * A helper method to calculate RMSE.
   */
  Eigen::Matrix<double, 4, 1, 0> CalculateRMSE(const std::vector<Eigen::Matrix<double, 4, 1, 0>> &estimations, 
                                            const std::vector<Eigen::Matrix<double, 4, 1, 0>> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

};

#endif  // TOOLS_H_
