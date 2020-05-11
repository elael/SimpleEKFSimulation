#pragma once

#include "Eigen/Dense"
#include "state.hpp"

template <int Y, int X>
class EKFUpdater {
 public:
  /**
   * Constructor
   * @param Inovation Inovation measurement function, usually Inovation(z,x) = z - h(x), where h( ) is the measurement function
   * @param J Jacobian of the measurement function
   * @param R Measurement covariance matrix function
   */
  EKFUpdater(InovationFunction<Y,X> Inovation, JacobianMap<Y,X> J, Eigen::Matrix<double, Y,Y> R);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   * @param x The state at k
   */
  State<X> Update(const Eigen::Matrix<double,Y,1> &z, const State<X>& x);

  private:
  // Inovation measurement function
  InovationFunction<Y, X> Inovation_;

  // Jacobian of the measurement function
  JacobianMap<Y, X> J_;

  // measurement covariance matrix function
  Eigen::Matrix<double, Y, Y> R_;
};

template <int Y, int X>
EKFUpdater<Y,X>::EKFUpdater(InovationFunction<Y,X> Inovation, JacobianMap<Y,X> J, Eigen::Matrix<double, Y,Y> R): Inovation_(Inovation), J_(J), R_(R)
{  
}

template <int Y, int X>
State<X> EKFUpdater<Y,X>::Update(const Eigen::Matrix<double,Y,1> &z, const State<X>& x){
  const auto inovation = Inovation_(z, x.mean);
  const auto H = J_(x.mean);
  const auto S = H*x.cov*H.transpose() + R_;
  const auto K = x.cov*H.transpose()*S.inverse();
  return State<X>{
    x.mean + K*inovation,
    (Eigen::Matrix<double, X, X>::Identity() - K*H)*x.cov
  };
}
