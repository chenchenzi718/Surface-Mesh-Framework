#ifndef SPECTRALQEM
#define SPECTRALQEM

#include "MeshViewerWidget.h"
#include<queue>
#include<map>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<algorithm>

#include<OpenMesh/Core/Utils/PropertyManager.hh>


class SpectralQEM {
public:
	SpectralQEM();
	~SpectralQEM();

	// 执行 arap 算法的示例代码
	void ARAP_kernel(Mesh& mesh);

	// 执行将 openmesh 转化为 eigen 数据供 igl 使用
	static void change_openmesh_to_eigen(Mesh mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	// 将 igl 数据转化为 openmesh 数据
	static void change_eigen_to_openmesh(Mesh& mesh, Eigen::MatrixXd& U);



public:
	// 设置一个结构体用来记录边的合并代价以及最优点
	struct EdgeCollapseInfo {
		Mesh::EdgeHandle eh;
		double cost;
		OpenMesh::Vec3d opt_point;
		int count; // 有些边进去优先队列时是更新前的信息，不作数
	};
	// 设置比较函数
	struct EdgeCompare {
		bool operator()(const EdgeCollapseInfo& a, const EdgeCollapseInfo& b) {
			// 根据MyClass对象中的value进行比较
			// 为了创建最小堆，这里应当返回true当第一个参数大于第二个参数
			return a.cost > b.cost;
		}
	};
	void Clear(std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare>& edge_que)
	{
		std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare> empty;
		swap(empty, edge_que);
	};

protected:



private:
	Eigen::MatrixXd V, U, bc;    // V 记录原始点的坐标，U 记录变换后坐标，bc 记录某些变换固定点坐标
	Eigen::MatrixXi F, b;        // F 记录面点索引，b 记录约束顶点索引

};


#endif // !SPECTRALQEM
