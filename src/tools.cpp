#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}


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

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  const double& vx = x_state(2);
  const double& vy = x_state(3);

  // squared distance
  float r2 = px*px + py*py;
  
  // If too close use velocity as direction
  if(r2 < 1e-10){
    px = vx;
    py = vy;
    r2 = px*px + py*py;
  }

  // If static return zero Jacobian
  if(r2 < 1e-10)    
      return MatrixXd::Zero(3, 4);  
  
  
  const float r = std::sqrt(r2);
  const float rho_px = px/r;
  const float rho_py = py/r;
  const float phi_px = - py/r2;
  const float phi_py = px / r2;
  
  const float proj = (vx*py - vy*px)/(r2*r);
  const float rhodot_px = py*proj;
  const float rhodot_py = (-px)*proj;
  
  // compute the Jacobian matrix
  Hj << rho_px, rho_py, 0, 0,
        phi_px, phi_py, 0, 0,
        rhodot_px, rhodot_py, rho_px, rho_py;

  return Hj;
}
