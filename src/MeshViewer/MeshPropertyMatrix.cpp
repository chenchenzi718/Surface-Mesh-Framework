#include "MeshPropertyMatrix.h"


double cotangent(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2){
	double cos_theta = v1.dot(v2);
	double sin_theta = v1.cross(v2).norm();
	return cos_theta / sin_theta;

};


void squared_edge_lengths(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& l2)
{
	l2.resize(F.rows(), 3);
	for (int i = 0; i < F.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			int idx0 = F(i, j);
			int idx1 = F(i, (j + 1) % 3);
			Eigen::Vector3d v0 = V.row(idx0);
			Eigen::Vector3d v1 = V.row(idx1);
			l2(i, j) = (v1 - v0).squaredNorm();
		}
	}
};

void squared_edge_lengths(const std::vector<Eigen::Vector3d>& V, Eigen::Vector3d& l2)
{
	for (int j = 0; j < 3; ++j) {
		auto v0 = V[j];
		auto v1 = V[(j + 1) % 3];
		l2(j) = (v1 - v0).squaredNorm();
	}
};


void doublearea(const Eigen::MatrixXd& l, Eigen::VectorXd& dblA) {
	dblA.resize(l.rows());
	for (int i = 0; i < l.rows(); ++i) {
		double a = l(i, 0);
		double b = l(i, 1);
		double c = l(i, 2);
		dblA(i) = 0.5 * std::sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)));
	}
};


double doublearea(const Eigen::Vector3d& l) {
	double a = l(0);
	double b = l(1);
	double c = l(2);
	return 0.5 * std::sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)));	
};


void computeCotangentWeights(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& C)
{
	// 计算边长的平方
	Eigen::MatrixXd l2;
	squared_edge_lengths(V, F, l2);

	// 计算边长
	Eigen::MatrixXd l = l2.array().sqrt();

	// 计算双面积
	Eigen::VectorXd dblA;
	doublearea(l, dblA);

	// 计算余弦权重
	// cosC = a^2+b^2-c^2/2ab, sinC=2S/ab, cotC = a^2+b^2-c^2/4S
	C.resize(F.rows(), 3);
	for (int i = 0; i < F.rows(); ++i) {
		C(i, 0) = (l2(i, 1) + l2(i, 2) - l2(i, 0)) / dblA(i) / 4.0;
		C(i, 1) = (l2(i, 2) + l2(i, 0) - l2(i, 1)) / dblA(i) / 4.0;
		C(i, 2) = (l2(i, 0) + l2(i, 1) - l2(i, 2)) / dblA(i) / 4.0;
	}
};

void computeCotangentWeights(const std::vector<Eigen::Vector3d>& V, const Eigen::Vector3d& l2, double dblA, Eigen::Vector3d& C)
{
	// 计算余弦权重
	// cosC = a^2+b^2-c^2/2ab, sinC=2S/ab, cotC = a^2+b^2-c^2/4S
	C(0) = (l2(1) + l2(2) - l2(0)) / dblA / 4.0;
	C(1) = (l2(2) + l2(0) - l2(1)) / dblA / 4.0;
	C(2) = (l2(0) + l2(1) - l2(2)) / dblA / 4.0;
};


void computeCotangentLaplacian(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& L)
{
	int n = V.rows();
	L.resize(n, n);

	// 计算余弦权重
	Eigen::MatrixXd C;
	computeCotangentWeights(V, F, C);

	std::vector<Eigen::Triplet<double>> triplets;

	for (int i = 0; i < F.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			int idx0 = F(i, j);
			int idx1 = F(i, (j + 1) % 3);

			double cot_ij = C(i, j);

			triplets.push_back(Eigen::Triplet<double>(idx0, idx1, cot_ij));
			triplets.push_back(Eigen::Triplet<double>(idx1, idx0, cot_ij));
			triplets.push_back(Eigen::Triplet<double>(idx0, idx0, -cot_ij));
			triplets.push_back(Eigen::Triplet<double>(idx1, idx1, -cot_ij));
		}
	}

	L.setFromTriplets(triplets.begin(), triplets.end());
}


void computeDiagonalMassMatrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& M) {
	int n = V.rows();
	Eigen::VectorXd dblA;
	Eigen::MatrixXd l2;

	// 计算边长的平方和双面积
	squared_edge_lengths(V, F, l2);

	// 计算边长
	Eigen::MatrixXd l = l2.array().sqrt();

	doublearea(l, dblA);

	// 计算每个顶点的质量
	Eigen::VectorXd vertex_mass = Eigen::VectorXd::Zero(n);
	for (int i = 0; i < F.rows(); ++i) {
		double area = dblA(i) / 2.0;
		for (int j = 0; j < 3; ++j) {
			vertex_mass(F(i, j)) += area / 3.0;
		}
	}

	// 构建对角质量矩阵
	std::vector<Eigen::Triplet<double>> triplets;
	for (int i = 0; i < n; ++i) {
		triplets.push_back(Eigen::Triplet<double>(i, i, vertex_mass(i)));
	}
	M.resize(n, n);
	M.setFromTriplets(triplets.begin(), triplets.end());
}