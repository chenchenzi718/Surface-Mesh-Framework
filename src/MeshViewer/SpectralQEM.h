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

	// ִ�� arap �㷨��ʾ������
	void ARAP_kernel(Mesh& mesh);

	// ִ�н� openmesh ת��Ϊ eigen ���ݹ� igl ʹ��
	static void change_openmesh_to_eigen(Mesh mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	// �� igl ����ת��Ϊ openmesh ����
	static void change_eigen_to_openmesh(Mesh& mesh, Eigen::MatrixXd& U);



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

protected:



private:
	Eigen::MatrixXd V, U, bc;    // V ��¼ԭʼ������꣬U ��¼�任�����꣬bc ��¼ĳЩ�任�̶�������
	Eigen::MatrixXi F, b;        // F ��¼���������b ��¼Լ����������

};


#endif // !SPECTRALQEM
