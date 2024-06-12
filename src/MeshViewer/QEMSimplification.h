#ifndef QEM
#define QEM

#include "MeshViewerWidget.h"
#include<queue>
#include<map>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<algorithm>


class QEMSimplification {
public:
	QEMSimplification(float ratio = 0.5f);
	~QEMSimplification();


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

	// r 为一个 整数
	void SetRatio(int r) {
		float ratio = (static_cast<float>(r)) / 100.f;
		this->ratio = std::max(0.f, std::min(1.f, ratio));
	};


	// 执行 QEM 算法
	void QEMSimplifyMesh(Mesh& mesh);


protected:
	// 设置每个 vertices 的 Q 矩阵
	void Cal_Q_ForVertex(const Mesh& mesh, Mesh::VertexHandle vertex_h);

	// 计算一条边合并之后的代价，并计算出一个最优点
	void Cal_Cost_ForEdge(const Mesh& mesh, Mesh::EdgeHandle eh);

	// 更新新的点处 Q 值以及临近边的重新计算
	void Update_Local_Variable(Mesh& mesh, Mesh::VertexHandle vertex_h);

	// 合并边
	bool EdgeCollapse(Mesh& mesh, const EdgeCollapseInfo& edge_collapse_info);


	float ratio;   // 简化等级

	double penalty;     // 设置边界边惩罚
	double cost_min;  // 如果最小的值弹出比这还大就不执行
	uint simplify_lowest_vertex_num;

	std::map<Mesh::VertexHandle, Eigen::Matrix4d> Q2v;
	std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare> edge_que;
	std::map<Mesh::EdgeHandle, int>each_edge_update_num_count;
};

#endif // !QEM
