#include "process_predictor.hpp"
#include "state.hpp"

// template <int X>
// ProcessPredictor<X>::ProcessPredictor(MatrixdFunction<X,X> F, MatrixdFunction<X,X>  Q): F_{F}, Q_{Q} {};

// /**
//  * Prediction Predicts the state and the state covariance
//  * using the process model
//  * @param delta_T Time between k and k+1 in s
//  */
// template <int X>
// State<X> ProcessPredictor<X>::Predict(long long delta_T, const State<X> &x){
//   auto F = F_(delta_T);
//   auto Q = Q_(delta_T);

//   return State{
//     F * x.mean,
//     F * x.cov * F.transpose() + Q
//   };
// }
