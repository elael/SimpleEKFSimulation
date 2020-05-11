#include "ekf_updater.hpp"
#include "Eigen/Dense"


// template <int Y, int X>
// EKFUpdater<Y,X>::EKFUpdater(InovationFunction<Y,X> Inovation, JacobianMap<Y,X> J, Eigen::Matrix<double, Y,Y> R): Inovation_(Inovation), J_(J), R_(R)
// {  
// }

// template <int Y, int X>
// State<X> EKFUpdater<Y,X>::Update(const Eigen::Matrix<double,Y,1> &z, const State<X>& x){
//   const auto inovation = Inovation_(z, x.mean);
//   const auto H = J_(x.mean);
//   const auto S = H*x.cov*H.transpose() + R_(x.mean);
//   const auto K = x.cov*H.transpose()*S.inverse();
//   return State{
//     x.mean + K*inovation,
//     (Eigen::Matrix<double, X, X>::Identity() - K*H)*x.cov
//   };
// }
