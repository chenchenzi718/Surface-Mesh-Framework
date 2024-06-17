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

	// 执行 arap 算法的示例代码
	void ARAP_kernel(Mesh& mesh);

	// 执行将 openmesh 转化为 eigen 数据供 igl 使用
	static void change_openmesh_to_eigen(Mesh mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	// 将 igl 数据转化为 openmesh 数据
	static void change_eigen_to_openmesh(Mesh& mesh, Eigen::MatrixXd& U);

    // 测试我们自己的 laplacian 矩阵与 igl 计算的
	void TestLaplacianIGLAndOurs(const Mesh& mesh);
	void TestMassMatrixIGLAndOurs(const Mesh& mesh);

    // SpetralQEM 主函数
    void SpectralQEMKernel(Mesh& mesh, size_t K = 100);


public:
    //  建立一个双向的 mapper，可以检查此时剩余的 vertex 内第 i 个对应的 vertex_id，也可以查询
    //  某个 vertex_id 排第几位
    class VertexMapper {
    public:
        // 添加顶点 ID
        void addVertex(int vertex_id) {
            vec_vid.push_back(vertex_id);
            position_in_vec_vid[vertex_id] = vec_vid.size() - 1;
        }

        // 通过索引获取顶点 ID
        int getVertexId(size_t index) const {
            if (index < vec_vid.size()) {
                return vec_vid[index];
            }
            else {
                throw std::out_of_range("Index out of range");
            }
        }

        // 通过顶点 ID 获取索引
        int getIndex(int vertex_id) const {
            auto it = position_in_vec_vid.find(vertex_id);
            if (it != position_in_vec_vid.end()) {
                return it->second;
            }
            else {
                throw std::invalid_argument("Vertex ID not found");
            }
        }

        // 删除顶点 ID
        void removeVertex(int vertex_id) {
            auto it = position_in_vec_vid.find(vertex_id);
            if (it != position_in_vec_vid.end()) {
                int index = it->second;
                vec_vid.erase(vec_vid.begin() + index);
                position_in_vec_vid.erase(it);

                // 更新 position_in_vec_vid 中的索引
                for (size_t i = index; i < vec_vid.size(); ++i) {
                    position_in_vec_vid[vec_vid[i]] = i;
                }
            }
            else {
                std::cerr << "Vertex ID not found in the vector" << std::endl;
            }
        }

        // 打印当前的顶点映射
        void printMapping() const {
            std::cout << "Current mapping: " << std::endl;
            for (size_t i = 0; i < vec_vid.size(); ++i) {
                std::cout << "Index " << i << " -> Vertex ID " << vec_vid[i] << std::endl;
            }
        }

        std::vector<int> vec_vid;    // 记录现存于 mesh 内的 vertex_id 的 vector，由小到大排列
        std::unordered_map<int, size_t> position_in_vec_vid;  // 根据 vertex id 得到其在 vec_vid 内的下标
    };


    // 设置一个结构体用来记录边的合并代价以及最优点
    struct EdgeCollapseInfo {
        Mesh::EdgeHandle eh;
        double cost;
        OpenMesh::Vec3d opt_point;
        int count; // 有些边进去优先队列时是更新前的信息，不作数

        SpMat L_bar; // 记录如果把这条 edge 进行 collapse 那么 L_bar 是什么样子的
        SpMat M_bar; // 记录如果把这条 edge 进行 collapse 那么 M_bar 是什么样子的
        Eigen::MatrixXd P_bar;
        VertexMapper mapper_bar;           // 记录如果把这条 edge 进行 collapse 那么 vertexmapper 是什么样子的
        uint collapse_signal;              // 这个 signal 反映的是该进行 heh 的 collapse 还是 heh_other 的 collapse
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
	// 将 mesh 内现存的 vertex id 在 vertex_ids 内由小到大排列，i1 是被删除的 vertex id
	void VertexReArrange(std::vector<int>& vertex_ids, int deleted_vertex_id);

    // 设置每个 vertices 的 Ev 矩阵
    void Cal_E_ForVertex(const Mesh& mesh, Eigen::MatrixXd P, Eigen::MatrixXd Z,
        SpMat M_bar, SpMat L_bar, Eigen::MatrixXd F,
        const VertexMapper& mapper_bar, std::map<int, double>& E2v);
    double Cal_E_ForVertex(const Mesh& mesh, Eigen::MatrixXd P, Eigen::MatrixXd Z,
        SpMat M_bar, SpMat L_bar, Eigen::MatrixXd F, const VertexMapper& mapper_bar, int index);

    // 计算每个 eh 上的能量
    void Cal_Cost_ForEdge(Mesh mesh,
        const SpMat& L, const SpMat& M, const Eigen::MatrixXd& P, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F,
        const VertexMapper& mapper,
        const Mesh::EdgeHandle eh);

    // 寻找一条边的 1 邻域的点集合
    void find_edge_1_ring(Mesh mesh, Mesh::EdgeHandle eh, std::vector<Mesh::EdgeHandle>& one_ring_pairs);

    // 根据 mapper, mapper_bar, mesh_bar, L, M 去构建新的 L_bar, M_bar
    void construct_new_LM(Mesh mesh, Mesh mesh_bar, const VertexMapper& mapper, const VertexMapper& mapper_bar,
        const SpMat& L, const SpMat& M, const Eigen::MatrixXd& P, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F,
        const std::vector<Mesh::EdgeHandle>& one_ring_pairs,
        int vidx_remain, int vidx_delete, EdgeCollapseInfo& info);

    // 建立 restriction matrix Q
    void construct_Q(const VertexMapper& mapper, const VertexMapper& mapper_bar, 
        int vidx_remain, int vidx_delete, SpMat& Q);

    void EdgeCollapse(Mesh& mesh, const EdgeCollapseInfo& edge_collapse_info, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F);
    void Update_Local_Variable(Mesh& mesh, 
        const SpMat& L, const SpMat& M, const Eigen::MatrixXd& P, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F,
        const VertexMapper& mapper,
        Mesh::VertexHandle vh);


private:
	Eigen::MatrixXd V, U, bc;    // V 记录原始点的坐标，U 记录变换后坐标，bc 记录某些变换固定点坐标
	Eigen::MatrixXi F, b;        // F 记录面点索引，b 记录约束顶点索引

    std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare> edge_que;
    std::map<Mesh::EdgeHandle, int>each_edge_update_num_count;

};


#endif // !SPECTRALQEM
