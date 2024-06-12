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
	// ����һ���ṹ��������¼�ߵĺϲ������Լ����ŵ�
	struct EdgeCollapseInfo {
		Mesh::EdgeHandle eh;
		double cost;
		OpenMesh::Vec3d opt_point;
		int count; // ��Щ�߽�ȥ���ȶ���ʱ�Ǹ���ǰ����Ϣ��������
	};
	// ���ñȽϺ���
	struct EdgeCompare {
		bool operator()(const EdgeCollapseInfo& a, const EdgeCollapseInfo& b) {
			// ����MyClass�����е�value���бȽ�
			// Ϊ�˴�����С�ѣ�����Ӧ������true����һ���������ڵڶ�������
			return a.cost > b.cost;
		}
	};
	void Clear(std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare>& edge_que)
	{
		std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare> empty;
		swap(empty, edge_que);
	};

	// r Ϊһ�� ����
	void SetRatio(int r) {
		float ratio = (static_cast<float>(r)) / 100.f;
		this->ratio = std::max(0.f, std::min(1.f, ratio));
	};


	// ִ�� QEM �㷨
	void QEMSimplifyMesh(Mesh& mesh);


protected:
	// ����ÿ�� vertices �� Q ����
	void Cal_Q_ForVertex(const Mesh& mesh, Mesh::VertexHandle vertex_h);

	// ����һ���ߺϲ�֮��Ĵ��ۣ��������һ�����ŵ�
	void Cal_Cost_ForEdge(const Mesh& mesh, Mesh::EdgeHandle eh);

	// �����µĵ㴦 Q ֵ�Լ��ٽ��ߵ����¼���
	void Update_Local_Variable(Mesh& mesh, Mesh::VertexHandle vertex_h);

	// �ϲ���
	bool EdgeCollapse(Mesh& mesh, const EdgeCollapseInfo& edge_collapse_info);


	float ratio;   // �򻯵ȼ�

	double penalty;     // ���ñ߽�߳ͷ�
	double cost_min;  // �����С��ֵ�������⻹��Ͳ�ִ��
	uint simplify_lowest_vertex_num;

	std::map<Mesh::VertexHandle, Eigen::Matrix4d> Q2v;
	std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare> edge_que;
	std::map<Mesh::EdgeHandle, int>each_edge_update_num_count;
};

#endif // !QEM
