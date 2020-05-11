#include "Eigen/Dense"

/* 
 * Laser Configuration
 */

// Measurement matrix
static const auto H_laser_matrix = [](){
  Eigen::Matrix<double, 2, 4> tmp;
  tmp << 1, 0, 0, 0,
         0, 1, 0, 0;
  return tmp;
}();

// Inovation function from state and measurement
static const auto Inovation_laser = [](const Eigen::Vector2d& z, const Eigen::Vector4d& x){
  return z - H_laser_matrix*x;
};

// Jacobian fuction of the measurement function
static const auto J_laser = [](auto){
  return H_laser_matrix;
};

// Covariance matrix of the measurement noise
static const auto R_laser = [](){
  Eigen::Matrix2d tmp;
  tmp << 0.0225, 0,
         0, 0.0225;
   return tmp;
}();

/* 
 * Radar Configuration
 */

void normalize_rad(double& theta){
  constexpr double pi = 3.14159265358979323846;
  if (std::abs(theta) <= pi) return;
  if (theta > pi)
    theta -= pi;
  else
    theta += pi;

  normalize_rad(theta);
}

// Inovation function from state and measurement
static const auto Inovation_radar = [](const Eigen::Vector3d& z, const Eigen::Vector4d& x){
  Eigen::Vector3d inovation;
  const double& px = x(0);
  const double& py = x(1);
  const double& vx = x(2);
  const double& vy = x(3);
  const double radius = std::sqrt(px*px + py*py);
  inovation(0) = z(0) - radius;
  inovation(1) = z(1) - std::atan2(py, px);
  normalize_rad(inovation(1));
  inovation(2) = z(2) - (px*vx + py*vy)/radius;
  return inovation;
};

// Jacobian fuction of the measurement function
static const auto J_radar = [](const Eigen::Vector4d& x_state) -> Eigen::Matrix<double, 3, 4> {

  Eigen::Matrix<double, 3, 4> Hj;
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  const double& vx = x_state(2);
  const double& vy = x_state(3);

  // squared distance
  double r2 = px*px + py*py;
  
  // If too close use velocity as direction
  if(r2 < 1e-20){
    px = vx;
    py = vy;
    r2 = px*px + py*py;
  }

  // If static return zero Jacobian
  if(r2 < 1e-20)    
      return Eigen::Matrix<double, 3, 4>::Zero();  
  
  
  const double r = std::sqrt(r2);
  const double rho_px = px/r;
  const double rho_py = py/r;
  const double phi_px = - py/r2;
  const double phi_py = px / r2;
  
  const double proj = (vx*py - vy*px)/(r2*r);
  const double rhodot_px = py*proj;
  const double rhodot_py = (-px)*proj;
  
  // compute the Jacobian matrix
  Hj << rho_px, rho_py, 0, 0,
        phi_px, phi_py, 0, 0,
        rhodot_px, rhodot_py, rho_px, rho_py;

  return Hj;
};

// Covariance matrix of the measurement noise
static const auto R_radar = [](){
  Eigen::Matrix3d tmp;
  tmp << 0.09, 0, 0,
         0, 0.0009, 0,
         0, 0, 0.09;
   return tmp;
}();


/* 
 * Process Configuration
 */

// Transition Matrix function
static const auto F = [](double dt){
  static Eigen::Matrix4d tmp = Eigen::Matrix4d::Identity();
  tmp(0,2) = dt; 
  tmp(1,3) = dt;
  return tmp;
};

// Process Covariance function
static const auto Q = [](double dt){
  static const double noise_ax = 9;
  static const double noise_ay = 9;
  Eigen::Matrix4d tmp;
  float dt2 = dt*dt;
  float dt3 = dt2 * dt /2.0;
  float dt4 = dt3 * dt /2.0;
  tmp << dt4*noise_ax, 0, dt3*noise_ax, 0,
          0, dt4*noise_ay, 0, dt3*noise_ay,
          dt3*noise_ax, 0, dt2 * noise_ax, 0,
          0, dt3*noise_ay, 0, dt2 * noise_ay;
  return tmp;
};

