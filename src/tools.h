#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

namespace Tools {
  /**
   * A helper method to calculate RMSE.
   */
  Eigen::Matrix<double, 4, 1, 0> CalculateRMSE(const std::vector<Eigen::Matrix<double, 4, 1, 0>> &estimations, 
                                            const std::vector<Eigen::Matrix<double, 4, 1, 0>> &ground_truth);

};

#endif  // TOOLS_H_
