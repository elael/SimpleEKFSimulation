#pragma once

#include "Eigen/Dense"

template <int N, int M>
using MatrixdFunction = std::function<const Eigen::Matrix<double, N, M> (double)>;

template <int N, int M>
using InovationFunction = std::function<const Eigen::Matrix<double, N, 1> (const Eigen::Matrix<double, N, 1>&, const Eigen::Matrix<double, M, 1>&)>;

template <int N, int M>
using JacobianMap = std::function<const Eigen::Matrix<double, N, M> (const Eigen::Matrix<double, M, 1>&)>;

template <int N, int M>
using VectorMap = std::function<Eigen::Matrix<double, N, 1> (const Eigen::Matrix<double, M, 1>&)>;

// Simple N dimensional state
template<int N>
struct State
{
  Eigen::Matrix<double, N, 1> mean;
  Eigen::Matrix<double, N, N> cov;
};
