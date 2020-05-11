#pragma once

#include "Eigen/Dense"
#include <functional>
#include "state.hpp"

template <int X>
class ProcessPredictor {
 public:
  /**
   * Constructor
   * @param F Transition matrix function
   * @param Q Process covariance matrix function
   */
  ProcessPredictor(MatrixdFunction<X,X> F, MatrixdFunction<X,X>  Q);
  
  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   * @param x The state at k
   */
  State<X> Predict(double delta_T, const State<X> &x);

  // state transition matrix function
  MatrixdFunction<X,X> F_;

  // process covariance matrix function
  MatrixdFunction<X,X> Q_;

};

template <int X>
ProcessPredictor<X>::ProcessPredictor(MatrixdFunction<X,X> F, MatrixdFunction<X,X>  Q): F_{F}, Q_{Q} {};

template <int X>
State<X> ProcessPredictor<X>::Predict(double delta_T, const State<X> &x){
  auto F = F_(delta_T);
  auto Q = Q_(delta_T);

  return State<X>{
    F * x.mean,
    F * x.cov * F.transpose() + Q
  };
}
