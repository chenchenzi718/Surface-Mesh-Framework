#ifndef MESHPROP
#define MESHPROP

#include "MeshViewerWidget.h"
#include<Eigen/Dense>
#include<Eigen/Core>
#include <Eigen/Sparse>


// 计算两个向量之间的 cot 值，构造 laplace 矩阵
double cotangent(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);

// 计算了所有三角形每条边 2 范数的平方
void squared_edge_lengths(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& l2);
// 针对某个三角形的计算
void squared_edge_lengths(const std::vector<Eigen::Vector3d>& V, Eigen::Vector3d& l2);

// 计算了每个三角形两倍的面积
void doublearea(const Eigen::MatrixXd& l, Eigen::VectorXd& dblA);
double doublearea(const Eigen::Vector3d& l);

// 计算每个三角形的余切权重
void computeCotangentWeights(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& C);
void computeCotangentWeights(const std::vector<Eigen::Vector3d>& V, const Eigen::Vector3d& l2, double dblA, Eigen::Vector3d& C);

// 计算 mesh 的 laplacian 矩阵
void computeCotangentLaplacian(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& L);

// 计算 mesh 的 mass matrix
void computeDiagonalMassMatrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& M);


#endif // !MESHPROP
