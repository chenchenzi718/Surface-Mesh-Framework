#ifndef MESHPROP
#define MESHPROP

#include "MeshViewerWidget.h"
#include<Eigen/Dense>
#include<Eigen/Core>
#include <Eigen/Sparse>


// ������������֮��� cot ֵ������ laplace ����
double cotangent(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);

// ����������������ÿ���� 2 ������ƽ��
void squared_edge_lengths(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& l2);
// ���ĳ�������εļ���
void squared_edge_lengths(const std::vector<Eigen::Vector3d>& V, Eigen::Vector3d& l2);

// ������ÿ�����������������
void doublearea(const Eigen::MatrixXd& l, Eigen::VectorXd& dblA);
double doublearea(const Eigen::Vector3d& l);

// ����ÿ�������ε�����Ȩ��
void computeCotangentWeights(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& C);
void computeCotangentWeights(const std::vector<Eigen::Vector3d>& V, const Eigen::Vector3d& l2, double dblA, Eigen::Vector3d& C);

// ���� mesh �� laplacian ����
void computeCotangentLaplacian(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& L);

// ���� mesh �� mass matrix
void computeDiagonalMassMatrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& M);


#endif // !MESHPROP
