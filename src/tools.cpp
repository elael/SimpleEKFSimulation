#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Eigen::Matrix<double, 4, 1> Tools::CalculateRMSE(const vector<Eigen::Matrix<double, 4, 1, 0>> &estimations,
                              const vector<Eigen::Matrix<double, 4, 1, 0>> &ground_truth) {
  Eigen::Matrix<double, 4, 1> rmse = Eigen::Matrix<double, 4, 1>::Zero();
  // rmse << 0,0,0,0;

  // Early return for empty array
  if (estimations.size() == 0)
    return rmse;

  // estimations and ground_truth should always have the same size, otherwise ignore ground_truth extra values
  assert(estimations.size() == ground_truth.size());

  // accumulate squared residuals
  for (std::size_t i = 0; i < estimations.size(); ++i) {
      Eigen::Matrix<double, 4, 1> error = estimations[i] - ground_truth[i];
      rmse += static_cast<Eigen::Matrix<double, 4, 1>>(error.array()*error.array());
  }

  // calculate the mean
  rmse /= estimations.size();

  // calculate the squared root and return the result
  return rmse.array().sqrt();
}
