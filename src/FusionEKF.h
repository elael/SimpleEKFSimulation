#ifndef FusionEKF_H_
#define FusionEKF_H_

#include <fstream>
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "ekf_updater.hpp"
#include "process_predictor.hpp"
#include "measurement_package.h"
#include "tools.h"
#include "state.hpp"
#include "process_predictor.hpp"

class FusionEKF {
 public:
  /**
   * Constructor.
   */
  FusionEKF();

  /**
   * Run the whole flow of the Kalman Filter from here.
   */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
   * Kalman Filter update and prediction math lives in here.
   */
  State<4> state;

 private:
  /**
   * Kalman Filter update and prediction math lives in here.
   */
  ProcessPredictor<4> process;
  EKFUpdater<2, 4> laser_filter;
  EKFUpdater<3, 4> radar_filter;

  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  Eigen::Matrix2d R_laser_;
  Eigen::Matrix3d R_radar_;
  Eigen::Matrix<double, 2, 4> H_laser_;
  Eigen::Matrix<double, 3, 4> Hj_;
};

#endif // FusionEKF_H_
