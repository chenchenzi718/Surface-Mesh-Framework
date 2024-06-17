#ifndef SPECTRALQEM
#define SPECTRALQEM

#include "MeshViewerWidget.h"
#include<queue>
#include<map>
#include<unordered_map>
#include <unordered_set>

#include<Eigen/Dense>
#include<Eigen/Core>
#include <Eigen/Sparse>
#include<algorithm>

#include"MeshPropertyMatrix.h"

typedef Eigen::SparseMatrix<double> SpMat;


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

    // ���������Լ��� laplacian ������ igl �����
	void TestLaplacianIGLAndOurs(const Mesh& mesh);
	void TestMassMatrixIGLAndOurs(const Mesh& mesh);

    // SpetralQEM ������
    void SpectralQEMKernel(Mesh& mesh, size_t K = 100);


public:
    //  ����һ��˫��� mapper�����Լ���ʱʣ��� vertex �ڵ� i ����Ӧ�� vertex_id��Ҳ���Բ�ѯ
    //  ĳ�� vertex_id �ŵڼ�λ
    class VertexMapper {
    public:
        // ��Ӷ��� ID
        void addVertex(int vertex_id) {
            vec_vid.push_back(vertex_id);
            position_in_vec_vid[vertex_id] = vec_vid.size() - 1;
        }

        // ͨ��������ȡ���� ID
        int getVertexId(size_t index) const {
            if (index < vec_vid.size()) {
                return vec_vid[index];
            }
            else {
                throw std::out_of_range("Index out of range");
            }
        }

        // ͨ������ ID ��ȡ����
        int getIndex(int vertex_id) const {
            auto it = position_in_vec_vid.find(vertex_id);
            if (it != position_in_vec_vid.end()) {
                return it->second;
            }
            else {
                throw std::invalid_argument("Vertex ID not found");
            }
        }

        // ɾ������ ID
        void removeVertex(int vertex_id) {
            auto it = position_in_vec_vid.find(vertex_id);
            if (it != position_in_vec_vid.end()) {
                int index = it->second;
                vec_vid.erase(vec_vid.begin() + index);
                position_in_vec_vid.erase(it);

                // ���� position_in_vec_vid �е�����
                for (size_t i = index; i < vec_vid.size(); ++i) {
                    position_in_vec_vid[vec_vid[i]] = i;
                }
            }
            else {
                std::cerr << "Vertex ID not found in the vector" << std::endl;
            }
        }

        // ��ӡ��ǰ�Ķ���ӳ��
        void printMapping() const {
            std::cout << "Current mapping: " << std::endl;
            for (size_t i = 0; i < vec_vid.size(); ++i) {
                std::cout << "Index " << i << " -> Vertex ID " << vec_vid[i] << std::endl;
            }
        }

        std::vector<int> vec_vid;    // ��¼�ִ��� mesh �ڵ� vertex_id �� vector����С��������
        std::unordered_map<int, size_t> position_in_vec_vid;  // ���� vertex id �õ����� vec_vid �ڵ��±�
    };


    // ����һ���ṹ��������¼�ߵĺϲ������Լ����ŵ�
    struct EdgeCollapseInfo {
        Mesh::EdgeHandle eh;
        double cost;
        OpenMesh::Vec3d opt_point;
        int count; // ��Щ�߽�ȥ���ȶ���ʱ�Ǹ���ǰ����Ϣ��������

        SpMat L_bar; // ��¼��������� edge ���� collapse ��ô L_bar ��ʲô���ӵ�
        SpMat M_bar; // ��¼��������� edge ���� collapse ��ô M_bar ��ʲô���ӵ�
        Eigen::MatrixXd P_bar;
        VertexMapper mapper_bar;           // ��¼��������� edge ���� collapse ��ô vertexmapper ��ʲô���ӵ�
        uint collapse_signal;              // ��� signal ��ӳ���Ǹý��� heh �� collapse ���� heh_other �� collapse
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
	// �� mesh ���ִ�� vertex id �� vertex_ids ����С�������У�i1 �Ǳ�ɾ���� vertex id
	void VertexReArrange(std::vector<int>& vertex_ids, int deleted_vertex_id);

    // ����ÿ�� vertices �� Ev ����
    void Cal_E_ForVertex(const Mesh& mesh, Eigen::MatrixXd P, Eigen::MatrixXd Z,
        SpMat M_bar, SpMat L_bar, Eigen::MatrixXd F,
        const VertexMapper& mapper_bar, std::map<int, double>& E2v);
    double Cal_E_ForVertex(const Mesh& mesh, Eigen::MatrixXd P, Eigen::MatrixXd Z,
        SpMat M_bar, SpMat L_bar, Eigen::MatrixXd F, const VertexMapper& mapper_bar, int index);

    // ����ÿ�� eh �ϵ�����
    void Cal_Cost_ForEdge(Mesh mesh,
        const SpMat& L, const SpMat& M, const Eigen::MatrixXd& P, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F,
        const VertexMapper& mapper,
        const Mesh::EdgeHandle eh);

    // Ѱ��һ���ߵ� 1 ����ĵ㼯��
    void find_edge_1_ring(Mesh mesh, Mesh::EdgeHandle eh, std::vector<Mesh::EdgeHandle>& one_ring_pairs);

    // ���� mapper, mapper_bar, mesh_bar, L, M ȥ�����µ� L_bar, M_bar
    void construct_new_LM(Mesh mesh, Mesh mesh_bar, const VertexMapper& mapper, const VertexMapper& mapper_bar,
        const SpMat& L, const SpMat& M, const Eigen::MatrixXd& P, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F,
        const std::vector<Mesh::EdgeHandle>& one_ring_pairs,
        int vidx_remain, int vidx_delete, EdgeCollapseInfo& info);

    // ���� restriction matrix Q
    void construct_Q(const VertexMapper& mapper, const VertexMapper& mapper_bar, 
        int vidx_remain, int vidx_delete, SpMat& Q);

    void EdgeCollapse(Mesh& mesh, const EdgeCollapseInfo& edge_collapse_info, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F);
    void Update_Local_Variable(Mesh& mesh, 
        const SpMat& L, const SpMat& M, const Eigen::MatrixXd& P, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F,
        const VertexMapper& mapper,
        Mesh::VertexHandle vh);


private:
	Eigen::MatrixXd V, U, bc;    // V ��¼ԭʼ������꣬U ��¼�任�����꣬bc ��¼ĳЩ�任�̶�������
	Eigen::MatrixXi F, b;        // F ��¼���������b ��¼Լ����������

    std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare> edge_que;
    std::map<Mesh::EdgeHandle, int>each_edge_update_num_count;

};


#endif // !SPECTRALQEM
