#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include "fusion_ekf_configuration.hpp"

using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

FusionEKF::FusionEKF(): 
    process(F, Q), 
    laser_filter(Inovation_laser, J_laser, R_laser), 
    radar_filter(Inovation_radar, J_radar, R_radar),
    is_initialized_(false)
{
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {

    state.cov = Eigen::Matrix4d::Random().cwiseAbs() * 1e-3; // add some uncertanties
    state.cov *= state.cov.transpose(); // cov are symmetrical
    state.cov.block<2,2>(2,2) += Eigen::Matrix2d::Identity() * 1e5; // add velocity uncertanty

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      // and initialize state.
      const double& rho = measurement_pack.raw_measurements_(0);
      const double& theta = measurement_pack.raw_measurements_(1);
      state.mean << rho*std::cos(theta), rho*std::sin(theta), 0, 0;
      Eigen::Matrix2d J;
      J << std::cos(theta), - rho * std::sin(theta),
           std::sin(theta),   rho * std::cos(theta);

      state.cov.block<2,2>(0,0) += J * R_radar.block<2,2>(0,0) * J.transpose();
    }
    else{
      state.mean << measurement_pack.raw_measurements_, 0, 0;
      state.cov.block<2,2>(0,0) += R_laser;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /**
   * Prediction
   */

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  state = process.Predict(dt, state);

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
    state = radar_filter.Update(Eigen::Ref<const Eigen::Vector3d>(measurement_pack.raw_measurements_), state);
  else 
    state = laser_filter.Update(Eigen::Ref<const Eigen::Vector2d>(measurement_pack.raw_measurements_), state);
    
}
