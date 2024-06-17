#include"SpectralQEM.h"

#include <igl/per_vertex_normals.h>
#include <igl/arap.h>
#include <OpenMesh/Core/Utils/PropertyManager.hh>



SpectralQEM::SpectralQEM()
{

}


SpectralQEM::~SpectralQEM()
{

}


void SpectralQEM::ARAP_kernel(Mesh& mesh)
{
    change_openmesh_to_eigen(mesh, this->V, this->F);

    // 计算顶点法线
    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);

    // 设置初始顶点位置
    U = V;

    // 定义变形后的控制点
    b.resize(2, 1); // 假设我们有两个控制点
    b << 1, 3;      // 约束顶点索引

    Eigen::MatrixXd bc(2, 3);  // bc 的行数应与 b 的长度相同
    bc.row(0) = V.row(1) + Eigen::RowVector3d(0.5, 0, 0); // 将顶点 1 向右移动 0.5 个单位
    bc.row(1) = V.row(3) + Eigen::RowVector3d(0.3, 0, 0); // 将顶点 3 向右移动 0.5 个单位

    // 设置 ARAP 算法参数
    igl::ARAPData arap_data;
    arap_data.with_dynamics = false;
    arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;

    // 预计算 ARAP
    arap_data.max_iter = 10;
    igl::arap_precomputation(V, F, V.cols(), b, arap_data);

    // 迭代求解 ARAP 变形
    igl::arap_solve(bc, arap_data, U);

    change_eigen_to_openmesh(mesh, U);

    
}


void SpectralQEM::change_openmesh_to_eigen(Mesh mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    V.resize(mesh.n_vertices(), 3);
    F.resize(mesh.n_faces(), 3);

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
        const auto& point = mesh.point(*v_it);
        V(v_it->idx(), 0) = point[0];
        V(v_it->idx(), 1) = point[1];
        V(v_it->idx(), 2) = point[2];
    }

    for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
    {
        auto fv_it = mesh.cfv_iter(*f_it);
        F(f_it->idx(), 0) = fv_it->idx(); ++fv_it;
        F(f_it->idx(), 1) = fv_it->idx(); ++fv_it;
        F(f_it->idx(), 2) = fv_it->idx();
    }
}


void SpectralQEM::change_eigen_to_openmesh(Mesh& mesh, Eigen::MatrixXd& U)
{
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
        mesh.set_point(*v_it, Mesh::Point(U(v_it->idx(), 0), U(v_it->idx(), 1), U(v_it->idx(), 2)));
    }
}


void SpectralQEM::TestLaplacianIGLAndOurs(const Mesh& mesh)
{
    change_openmesh_to_eigen(mesh, this->V, this->F);

    SpMat L_igl;
    igl::cotmatrix(V, F, L_igl);

    // 计算手动的拉普拉斯矩阵
    SpMat L;
    computeCotangentLaplacian(this->V, this->F, L);

    auto compareMatrices = [](const SpMat& A, const SpMat& B, double threshold)
    {
        int num_differences = 0;
        for (int k = 0; k < A.outerSize(); ++k) {
            for (SpMat::InnerIterator it(A, k); it; ++it) {
                double diff = std::abs(it.value() - B.coeff(it.row(), it.col()));
                if (diff > threshold) {
                    ++num_differences;
                    std::cout << "Difference at (" << it.row() << ", " << it.col() << "): "
                        << it.value() << " vs " << B.coeff(it.row(), it.col()) << " (diff: " << diff << ")\n";
                }
            }
        }
        return num_differences;
    };

    // 比较两个拉普拉斯矩阵
    double threshold = 1e-6;
    int num_differences = compareMatrices(L, L_igl, threshold);

    std::cout << "Number of differing elements: " << num_differences << std::endl;
}


void SpectralQEM::TestMassMatrixIGLAndOurs(const Mesh& mesh)
{
    change_openmesh_to_eigen(mesh, this->V, this->F);

    // 计算手动的对角质量矩阵
    SpMat M;
    computeDiagonalMassMatrix(V, F, M);

    // 使用 libigl 计算对角质量矩阵
    SpMat M_igl;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M_igl);

    auto compareMatrices = [](const SpMat& A, const SpMat& B, double threshold)
    {
        int num_differences = 0;
        for (int k = 0; k < A.outerSize(); ++k) {
            for (SpMat::InnerIterator it(A, k); it; ++it) {
                double diff = std::abs(it.value() - B.coeff(it.row(), it.col()));
                if (diff > threshold) {
                    ++num_differences;
                    std::cout << "Difference at (" << it.row() << ", " << it.col() << "): "
                        << it.value() << " vs " << B.coeff(it.row(), it.col()) << " (diff: " << diff << ")\n";
                }
            }
        }
        return num_differences;
    };

    // 比较两个拉普拉斯矩阵
    double threshold = 1e-6;
    int num_differences = compareMatrices(M, M_igl, threshold);

    std::cout << "Number of differing elements: " << num_differences << std::endl;
}


void SpectralQEM::VertexReArrange(std::vector<int>& vertex_ids, int deleted_vertex_id)
{
    // 找到值为 value 的元素位置
    auto it = std::find(vertex_ids.begin(), vertex_ids.end(), deleted_vertex_id);
    if (it != vertex_ids.end()) {
        // 删除该元素
        vertex_ids.erase(it);
    }
    else {
        std::cerr << "Value not found in the vector" << std::endl;
    }
}


void SpectralQEM::SpectralQEMKernel(Mesh& mesh, size_t K)
{
    change_openmesh_to_eigen(mesh, this->V, this->F);
    uint mesh_vertex_num = mesh.n_vertices();
    uint simplify_mesh_vertex = static_cast<uint>(0.5 * mesh_vertex_num);
    if (simplify_mesh_vertex < 3)
        simplify_mesh_vertex = 3;

    // 手动计算的拉普拉斯矩阵与对角质量矩阵
    SpMat L, M;
    computeCotangentLaplacian(V, F, L);
    computeDiagonalMassMatrix(V, F, M);

    // 创建 VertexMapper 
    VertexMapper mapper;
    for (auto vh : mesh.vertices()) {
        mapper.addVertex(vh.idx());
    }

    size_t n = mesh.n_vertices();
    // 创建初始的 fine-to-coarse restriction matrix
    Eigen::MatrixXd P = Eigen::MatrixXd::Identity(n, n);


    // 创建一个 EigenSolver 对象
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(L);
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Eigen decomposition failed!" << std::endl;
        return;
    }
    // 获取前 K 个特征向量
    Eigen::MatrixXd F = eigensolver.eigenvectors().leftCols(K);


    // 计算 M 的逆矩阵
    SpMat M_inv(n, n);
    for (int k = 0; k < M.outerSize(); ++k) {
        for (SpMat::InnerIterator it(M, k); it; ++it) {
            M_inv.insert(it.row(), it.col()) = 1.0 / it.value();
        }
    }

    // 预计算矩阵 Z
    Eigen::MatrixXd Z = (M_inv * L) * F;
    SpMat M_bar = M;
    SpMat L_bar = L;

    // 从 vertex.idx() 指向它们的能量
    std::map<int, double> E2v;  
    // 初始每个点的能量，实际上都是 0
    for (int i = 0; i < n; ++i) {
        E2v[i] = 0.0;
    }

    // 计算每条边的 cost
    for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
    {
        each_edge_update_num_count[*e_it] = 0;
        Cal_Cost_ForEdge(mesh, L_bar, M_bar, P, Z, F, mapper, *e_it);
    }

    int collapse_num = 0;
    // 进行循环
    while (collapse_num < (mesh_vertex_num - simplify_mesh_vertex) && !edge_que.empty())
    {
        // 取出最小的 cost 边
        EdgeCollapseInfo edgeinfo = edge_que.top();
        edge_que.pop();

        auto eh = edgeinfo.eh;
        if (mesh.status(eh).deleted())
            continue;
        if (edgeinfo.count != each_edge_update_num_count[eh])
            continue;

        EdgeCollapse(mesh, edgeinfo, Z, F);
        collapse_num++;
    }

    Clear(edge_que); // 如果不把 handle 清理干净，garbage 收集之后就报错
    each_edge_update_num_count.clear();

    mesh.garbage_collection();
}


void SpectralQEM::EdgeCollapse(Mesh& mesh, const EdgeCollapseInfo& edge_collapse_info, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F)
{
    auto eh = edge_collapse_info.eh;
    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh);
    Mesh::HalfedgeHandle heh_other = mesh.opposite_halfedge_handle(heh);

    bool collapse_ok = false;

    auto from_p = mesh.from_vertex_handle(heh);
    auto to_p = mesh.to_vertex_handle(heh);

    auto sp = edge_collapse_info.opt_point;
    Mesh::VertexHandle v_result;


    if (!edge_collapse_info.collapse_signal)
    {
        mesh.set_point(to_p, sp);
        mesh.collapse(heh);
        v_result = to_p;
    }
    else
    {
        mesh.set_point(from_p, sp);
        mesh.collapse(heh_other);
        v_result = from_p;
    }

    Update_Local_Variable(mesh, edge_collapse_info.L_bar,edge_collapse_info.M_bar, edge_collapse_info.P_bar, Z, F, 
        edge_collapse_info.mapper_bar, v_result);
}


void SpectralQEM::Update_Local_Variable(Mesh& mesh, 
    const SpMat& L, const SpMat& M, const Eigen::MatrixXd& P, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F, 
    const VertexMapper& mapper,
    Mesh::VertexHandle vh)
{
    for (auto ve_it = mesh.ve_iter(vh); ve_it.is_valid(); ++ve_it)
    {
        each_edge_update_num_count[*ve_it]++;
        Cal_Cost_ForEdge(mesh, L, M, P, Z, F, mapper, *ve_it);
    }
}


void SpectralQEM::Cal_E_ForVertex(const Mesh& mesh, Eigen::MatrixXd P, Eigen::MatrixXd Z, 
    SpMat M_bar, SpMat L_bar, Eigen::MatrixXd F, 
    const VertexMapper &mapper_bar, std::map<int, double>& E2v)
{
    size_t n_bar = L_bar.rows();

    // 计算 M 的逆矩阵
    SpMat M_bar_inv(n_bar, n_bar);
    for (int k = 0; k < M_bar.outerSize(); ++k) {
        for (SpMat::InnerIterator it(M_bar, k); it; ++it) {
            M_bar_inv.insert(it.row(), it.col()) = 1.0 / it.value();
        }
    }
    M_bar_inv.makeCompressed();

    // 计算 \tilde{M}^{-1} \tilde{L} P F
    Eigen::MatrixXd M_inv_LPF = M_bar_inv * L_bar * P * F;

    //  计算每个顶点的 E2v
    for (size_t i = 0; i < n_bar; ++i)
    {
        double Mv = M_bar.coeff(i, i);
        Eigen::VectorXd row_diff = P.row(i) * Z - M_inv_LPF.row(i);
        E2v[mapper_bar.getVertexId(i)] = Mv * row_diff.squaredNorm();
    }
}


double SpectralQEM::Cal_E_ForVertex(const Mesh& mesh, Eigen::MatrixXd P, Eigen::MatrixXd Z,
    SpMat M_bar, SpMat L_bar, Eigen::MatrixXd F, const VertexMapper& mapper_bar, int v_index)
{
    size_t n_bar = L_bar.rows();

    // 计算 M 的逆矩阵
    SpMat M_bar_inv(n_bar, n_bar);
    for (int k = 0; k < M_bar.outerSize(); ++k) {
        for (SpMat::InnerIterator it(M_bar, k); it; ++it) {
            M_bar_inv.insert(it.row(), it.col()) = 1.0 / it.value();
        }
    }
    M_bar_inv.makeCompressed();

    // 计算 \tilde{M}^{-1} \tilde{L} P F
    Eigen::MatrixXd M_inv_LPF = M_bar_inv * L_bar * P * F;

    //  计算每个顶点的 E2v
    int index = mapper_bar.getIndex(v_index);
    double Mv = M_bar.coeff(index, index);
    Eigen::VectorXd row_diff = P.row(index) * Z - M_inv_LPF.row(index);
    return Mv * row_diff.squaredNorm(); 
}


void SpectralQEM::Cal_Cost_ForEdge(Mesh mesh, 
    const SpMat& L, const SpMat& M, const Eigen::MatrixXd& P, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F,
    const VertexMapper& mapper, 
    const Mesh::EdgeHandle eh)
{
    EdgeCollapseInfo info;
    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh);
    Mesh::HalfedgeHandle heh_other = mesh.opposite_halfedge_handle(heh);

    Mesh::VertexHandle from_v = mesh.from_vertex_handle(heh);
    Mesh::VertexHandle to_v = mesh.to_vertex_handle(heh);
    auto from_p = mesh.point(from_v); 
    auto to_p = mesh.point(to_v);
    OpenMesh::Vec3d opt_point = (from_p + to_p) / 2.0;

    info.eh = eh;
    info.opt_point = opt_point;
    info.count = each_edge_update_num_count[eh];

    // 下面计算这条边的 cost 以及将它 collapse 之后的 L_bar,M_bar,mapper_bar
    Mesh mesh_bar = mesh;
    int vidx_remain = to_v.idx();
    int vidx_delete = from_v.idx();
    std::vector<Mesh::EdgeHandle> one_ring_pairs;
    find_edge_1_ring(mesh, eh, one_ring_pairs);

    if (mesh.is_collapse_ok(heh))
    {
        mesh_bar.set_point(to_v, opt_point);
        mesh_bar.collapse(heh);
        info.collapse_signal = 0;
    }
    else if (mesh.is_collapse_ok(heh_other))
    {
        vidx_remain = from_v.idx();
        vidx_delete = to_v.idx();

        mesh_bar.set_point(from_v, opt_point);
        mesh_bar.collapse(heh_other);
        info.collapse_signal = 1;
    }
    else
    {
        // 否则不将这条边加入可以 collapse 的行列
        return;
    }

    VertexMapper mapper_bar = mapper;
    mapper_bar.removeVertex(vidx_delete);
    info.mapper_bar = mapper_bar;

    // 根据 mapper, mapper_bar, mesh_bar, L, M 去构建新的 L_bar, M_bar
    construct_new_LM(mesh, mesh_bar, mapper, mapper_bar, L, M, P, Z, F, one_ring_pairs, vidx_remain, vidx_delete, info);
    edge_que.push(info);
}


void SpectralQEM::find_edge_1_ring(Mesh mesh, Mesh::EdgeHandle eh, std::vector<Mesh::EdgeHandle>& one_ring_pairs)
{
    // 获取这条边的两个顶点
    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
    Mesh::VertexHandle to_v = mesh.to_vertex_handle(heh);
    Mesh::VertexHandle from_v = mesh.from_vertex_handle(heh);

    // 遍历与这条边的两个顶点相邻的所有半边，找到所有相邻的面
    std::unordered_set<Mesh::FaceHandle> adjacent_faces;
    for (auto voh_it = mesh.voh_begin(from_v); voh_it != mesh.voh_end(from_v); ++voh_it) {
        Mesh::FaceHandle fh = mesh.face_handle(*voh_it);
        if (fh.is_valid()) { 
            adjacent_faces.insert(fh);
        }
    }
    for (auto voh_it = mesh.voh_begin(to_v); voh_it != mesh.voh_end(to_v); ++voh_it) {
        Mesh::FaceHandle fh = mesh.face_handle(*voh_it);
        if (fh.is_valid()) {
            adjacent_faces.insert(fh);
        }
    }

    // 对于每个相邻的面，找到除这条边的两个顶点外的其他顶点构成的边
    for (const auto& fh : adjacent_faces) {
        std::vector<Mesh::VertexHandle> vertices;
        for (auto fe_it = mesh.fe_begin(fh); fe_it != mesh.fe_end(fh); ++fe_it) {
            Mesh::EdgeHandle eh = *fe_it;          
            Mesh::HalfedgeHandle heh0 = mesh.halfedge_handle(eh, 0);
            Mesh::VertexHandle to_vh = mesh.to_vertex_handle(heh0);
            Mesh::VertexHandle from_vh = mesh.from_vertex_handle(heh0);

            if ((to_vh != from_v) && (to_vh != to_v) && (from_vh != from_v) && (from_vh != to_v)) {
                one_ring_pairs.push_back(eh);
            }
        }
    }
}


void SpectralQEM::construct_new_LM(Mesh mesh, Mesh mesh_bar, const VertexMapper& mapper, const VertexMapper& mapper_bar,
    const SpMat& L, const SpMat& M, const Eigen::MatrixXd& P, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& F,
    const std::vector<Mesh::EdgeHandle>& one_ring_pairs,
    int vidx_remain, int vidx_delete, EdgeCollapseInfo& info)
{
    std::vector<Eigen::Triplet<double>> triplets_L_bar, triplets_M_bar;
    SpMat L_bar, M_bar;
    int n_bar = mapper_bar.vec_vid.size();
    L_bar.resize(n_bar, n_bar);
    M_bar.resize(n_bar, n_bar);

    std::unordered_set<int> one_ring_vertices;
    for (const auto& eh : one_ring_pairs)
    {
        Mesh::HalfedgeHandle heh0 = mesh_bar.halfedge_handle(eh, 0);
        Mesh::VertexHandle to_vh = mesh_bar.to_vertex_handle(heh0);
        Mesh::VertexHandle from_vh = mesh_bar.from_vertex_handle(heh0);
        one_ring_vertices.insert(to_vh.idx());
        one_ring_vertices.insert(from_vh.idx());
    }

    for (int k = 0; k < L.outerSize(); ++k) {
        for (SpMat::InnerIterator it(L, k); it; ++it) {
            //  这反映的是 mapper 内的 index
            int index_1 = it.row();
            int index_2 = it.col();

            auto vidx_1 = mapper.getVertexId(index_1);
            auto vidx_2 = mapper.getVertexId(index_2);

            // 首先 vidx_1, vidx_2 内不能有 vidx_remain，vidx_delete ，否则重新算
            // 其次 vidx_1 = vidx_2 时它们不能落在 one_ring_pairs 内，否则重新算
            // 最后 vidx_1 不等于 vidx_2 时它们作为一个 pair 不能落在 one_ring_pairs 内，否则重新算
            if ((vidx_1 != vidx_remain) && (vidx_1 != vidx_delete) && (vidx_2 != vidx_remain) && (vidx_2 != vidx_delete))
            {
                // vidx_1, vidx_2 内此时不存在 vidx_remain, vidx_delete
                if (vidx_1 == vidx_2)
                {
                    // 如果一样，不能落在 one_ring_vertices 内
                    if (one_ring_vertices.count(vidx_1) == 0)
                    {
                        // 此时的 L 可以重复使用
                        int index_1_bar = mapper_bar.getIndex(vidx_1);
                        triplets_L_bar.push_back(Eigen::Triplet<double>(index_1_bar, index_1_bar, it.value()));
                    }
                }
                else
                {
                    if (one_ring_vertices.count(vidx_1) == 0 || one_ring_vertices.count(vidx_2) == 0)
                    {
                        // 此时的 L 可以重复使用
                        int index_1_bar = mapper_bar.getIndex(vidx_1);
                        int index_2_bar = mapper_bar.getIndex(vidx_2);
                        triplets_L_bar.push_back(Eigen::Triplet<double>(index_1_bar, index_2_bar, it.value()));
                    }
                }
            }
        }
    }
    // 特殊的边，点进行重新计算
    // 首先收集 one_ring_vertices 会涉及到的面
    std::unordered_set<Mesh::FaceHandle> face_set;
    // 记录每个面的双面积
	auto double_area = OpenMesh::makeTemporaryProperty<OpenMesh::FaceHandle, double>(mesh_bar);
    // 记录每个面的三个点的 idx
    auto vertex_id_in_face = OpenMesh::makeTemporaryProperty<OpenMesh::FaceHandle, Eigen::Vector3i>(mesh_bar);
    // 记录每个面的三个 cot 值，与 vertex 对应
    auto cot_in_face = OpenMesh::makeTemporaryProperty<OpenMesh::FaceHandle, Eigen::Vector3d>(mesh_bar);

    for (const auto& vidx : one_ring_vertices)
    {
        Mesh::VertexHandle vh = mesh_bar.vertex_handle(vidx);
        for (auto vf_iter = mesh_bar.vf_begin(vh); vf_iter != mesh_bar.vf_end(vh); ++vf_iter)
        {
            if((*vf_iter).is_valid())
                face_set.insert(*vf_iter);
        }
    }
    for (const auto& fh : face_set)
    {
        std::vector<Eigen::Vector3d>points;
        int count = 0;
        for (auto fv_iter = mesh_bar.fv_begin(fh); fv_iter != mesh_bar.fv_end(fh); fv_iter++)
        {
            auto vh = *fv_iter;
            vertex_id_in_face[fh][count] = vh.idx(); count++;
            auto p = mesh_bar.point(vh);
            points.push_back(Eigen::Vector3d(p[0], p[1], p[2]));
        }
        // 计算边长的平方
        Eigen::Vector3d l2;
        squared_edge_lengths(points, l2);
        // 计算边长
        Eigen::Vector3d l = l2.array().sqrt();
        double_area[fh] = doublearea(l);
        Eigen::Vector3d C;
        computeCotangentWeights(points, l2, double_area[fh], C);
        cot_in_face[fh] = C;
    }

    // 首先更新每一条 one-ring 边的 L
    for (const auto& eh : one_ring_pairs)
    {
        Mesh::HalfedgeHandle heh0 = mesh_bar.halfedge_handle(eh, 0);
        Mesh::HalfedgeHandle heh1 = mesh_bar.halfedge_handle(eh, 1);
        Mesh::VertexHandle to_vh = mesh_bar.to_vertex_handle(heh0);
        Mesh::VertexHandle from_vh = mesh_bar.from_vertex_handle(heh0);
        
        Mesh::FaceHandle fh = mesh_bar.face_handle(heh0);
        if (fh.is_valid())
        {
            auto vertex_id_arr = vertex_id_in_face[fh];
            int count = 0;
            for (count = 0; count < 3; count++)
            {
                if (vertex_id_arr[count] == from_vh.idx()) break;
            }
            auto index_1 = mapper_bar.getIndex(from_vh.idx());
            auto index_2 = mapper_bar.getIndex(to_vh.idx());
            triplets_L_bar.push_back(Eigen::Triplet<double>(index_1, index_2, cot_in_face[fh][count]));
            triplets_L_bar.push_back(Eigen::Triplet<double>(index_2, index_1, cot_in_face[fh][count]));
        }

        fh = mesh_bar.face_handle(heh1);
        if (fh.is_valid())
        {
            auto vertex_id_arr = vertex_id_in_face[fh];
            int count = 0;
            for (count = 0; count < 3; count++)
            {
                if (vertex_id_arr[count] == to_vh.idx()) break;
            }
            auto index_1 = mapper_bar.getIndex(from_vh.idx());
            auto index_2 = mapper_bar.getIndex(to_vh.idx());
            triplets_L_bar.push_back(Eigen::Triplet<double>(index_1, index_2, cot_in_face[fh][count]));
            triplets_L_bar.push_back(Eigen::Triplet<double>(index_2, index_1, cot_in_face[fh][count]));
        }
    }
    // 更新每一个 one-ring 上面的点的 L
    for (const auto& vidx : one_ring_vertices)
    {
        Mesh::VertexHandle vh = mesh_bar.vertex_handle(vidx);
        for (auto voh_iter = mesh_bar.voh_begin(vh); voh_iter != mesh_bar.voh_end(vh); ++voh_iter)
        {
            auto heh = *voh_iter;
            Mesh::VertexHandle to_vh = mesh_bar.to_vertex_handle(heh);
            Mesh::VertexHandle from_vh = mesh_bar.from_vertex_handle(heh);

            Mesh::FaceHandle fh = mesh_bar.face_handle(heh);
            if (fh.is_valid())
            {
                auto vertex_id_arr = vertex_id_in_face[fh];
                int count = 0;
                for (count = 0; count < 3; count++)
                {
                    if (vertex_id_arr[count] == from_vh.idx()) break;
                }
                auto index_1 = mapper_bar.getIndex(from_vh.idx());
                triplets_L_bar.push_back(Eigen::Triplet<double>(index_1, index_1, -cot_in_face[fh][count]));
            }
        }
    }
    // 更新从 vidx_remain 出去的所有边
    auto vh_remain = mesh_bar.vertex_handle(vidx_remain);
    for (auto ve_iter = mesh_bar.ve_begin(vh_remain); ve_iter != mesh_bar.ve_end(vh_remain); ++ve_iter)
    {
        auto eh = *ve_iter;
        Mesh::HalfedgeHandle heh0 = mesh_bar.halfedge_handle(eh, 0);
        Mesh::HalfedgeHandle heh1 = mesh_bar.halfedge_handle(eh, 1);
        Mesh::VertexHandle to_vh = mesh_bar.to_vertex_handle(heh0);
        Mesh::VertexHandle from_vh = mesh_bar.from_vertex_handle(heh0);

        Mesh::FaceHandle fh = mesh_bar.face_handle(heh0);
        if (fh.is_valid())
        {
            auto vertex_id_arr = vertex_id_in_face[fh];
            int count = 0;
            for (count = 0; count < 3; count++)
            {
                if (vertex_id_arr[count] == from_vh.idx()) break;
            }
            auto index_1 = mapper_bar.getIndex(from_vh.idx());
            auto index_2 = mapper_bar.getIndex(to_vh.idx());
            triplets_L_bar.push_back(Eigen::Triplet<double>(index_1, index_2, cot_in_face[fh][count]));
            triplets_L_bar.push_back(Eigen::Triplet<double>(index_2, index_1, cot_in_face[fh][count]));
        }

        fh = mesh_bar.face_handle(heh1);
        if (fh.is_valid())
        {
            auto vertex_id_arr = vertex_id_in_face[fh];
            int count = 0;
            for (count = 0; count < 3; count++)
            {
                if (vertex_id_arr[count] == to_vh.idx()) break;
            }
            auto index_1 = mapper_bar.getIndex(from_vh.idx());
            auto index_2 = mapper_bar.getIndex(to_vh.idx());
            triplets_L_bar.push_back(Eigen::Triplet<double>(index_1, index_2, cot_in_face[fh][count]));
            triplets_L_bar.push_back(Eigen::Triplet<double>(index_2, index_1, cot_in_face[fh][count]));
        }
    }
    // 更新 vidx_remain 对应的值
    for (auto voh_iter = mesh_bar.voh_begin(vh_remain); voh_iter != mesh_bar.voh_end(vh_remain); ++voh_iter)
    {
        auto heh = *voh_iter;
        Mesh::VertexHandle to_vh = mesh_bar.to_vertex_handle(heh);
        Mesh::VertexHandle from_vh = mesh_bar.from_vertex_handle(heh);

        Mesh::FaceHandle fh = mesh_bar.face_handle(heh);
        if (fh.is_valid())
        {
            auto vertex_id_arr = vertex_id_in_face[fh];
            int count = 0;
            for (count = 0; count < 3; count++)
            {
                if (vertex_id_arr[count] == from_vh.idx()) break;
            }
            auto index_1 = mapper_bar.getIndex(from_vh.idx());
            triplets_L_bar.push_back(Eigen::Triplet<double>(index_1, index_1, -cot_in_face[fh][count]));
        }
    }
    L_bar.setFromTriplets(triplets_L_bar.begin(), triplets_L_bar.end());
    info.L_bar = L_bar;

    ///////////////
    for (int k = 0; k < M.outerSize(); ++k) {
        for (SpMat::InnerIterator it(M, k); it; ++it) {
            //  这反映的是 mapper 内的 index
            int index = it.row();
            auto vidx = mapper.getVertexId(index);

            // vidx 不能落在 one_ring_vertices,vidx_remain,vidx_delete  内，否则重新算
            if (vidx != vidx_remain && vidx != vidx_delete)
            {
                // 如果一样，不能落在 one_ring_vertices 内
                if (one_ring_vertices.count(vidx) == 0)
                {
                    // 此时的 L 可以重复使用
                    int index_bar = mapper_bar.getIndex(vidx);
                    triplets_M_bar.push_back(Eigen::Triplet<double>(index_bar, index_bar, it.value()));
                }             
            }
        }
    }
    // 更新 one_ring_vertices 上的
    for (const auto& vidx : one_ring_vertices)
    {
        Mesh::VertexHandle vh = mesh_bar.vertex_handle(vidx);
        double mass = 0.0;
        for (auto vf_iter = mesh_bar.vf_begin(vh); vf_iter != mesh_bar.vf_end(vh); ++vf_iter)
        {
            if ((*vf_iter).is_valid())
                mass += double_area[*vf_iter] / 6.0;
        }
        auto index_bar = mapper_bar.getIndex(vidx);
        triplets_M_bar.push_back(Eigen::Triplet<double>(index_bar, index_bar, mass));
    }
    // 更新 vidx_remain 的值
    double mass = 0.0;
    for (auto vf_iter = mesh_bar.vf_begin(vh_remain); vf_iter != mesh_bar.vf_end(vh_remain); ++vf_iter)
    {
        if ((*vf_iter).is_valid())
        {
            mass += double_area[*vf_iter] / 6.0;
        }
    }
    auto index_bar = mapper_bar.getIndex(vidx_remain);
    triplets_M_bar.push_back(Eigen::Triplet<double>(index_bar, index_bar, mass));
    
    M_bar.setFromTriplets(triplets_M_bar.begin(), triplets_M_bar.end());
    info.M_bar = M_bar;

    // 计算边的 cost
    SpMat Q;
    construct_Q(mapper, mapper_bar, vidx_remain, vidx_delete, Q);
    auto P_bar = Q * P;
    double cost_pre = 0.0;
    double cost_now = 0.0;
    for (const auto& vidx : one_ring_vertices)
    {
        cost_pre += Cal_E_ForVertex(mesh, P, Z, M, L, F, mapper, vidx);
        cost_now += Cal_E_ForVertex(mesh_bar, P_bar, Z, M_bar, L_bar, F, mapper_bar, vidx);
    }
    cost_pre += (Cal_E_ForVertex(mesh, P, Z, M, L, F, mapper, vidx_remain) + Cal_E_ForVertex(mesh, P, Z, M, L, F, mapper, vidx_delete));
    cost_now += Cal_E_ForVertex(mesh_bar, P_bar, Z, M_bar, L_bar, F, mapper_bar, vidx_remain);
    info.P_bar = P_bar;
    info.cost = cost_now - cost_pre;
}

void SpectralQEM::construct_Q(const VertexMapper& mapper, const VertexMapper& mapper_bar, 
    int vidx_remain, int vidx_delete,
    SpMat& Q)
{
    std::vector<Eigen::Triplet<double>> triplets;
    int n_pre = mapper.vec_vid.size();
    int n_now = mapper_bar.vec_vid.size();
    Q.resize(n_now, n_pre);

    for (const auto& pair : mapper_bar.position_in_vec_vid)
    {
        if (pair.first != vidx_remain)
        {
            int index_now = pair.second;
            int index_pre = mapper.getIndex(pair.first);
            triplets.push_back(Eigen::Triplet<double>(index_now, index_pre, 1.0));
        }
        else {
            int index_now = pair.second;
            int index_pre = mapper.getIndex(vidx_remain);
            int index_pre_2 = mapper.getIndex(vidx_delete);
            triplets.push_back(Eigen::Triplet<double>(index_now, index_pre, 0.5));
            triplets.push_back(Eigen::Triplet<double>(index_now, index_pre_2, 0.5));
        }
    }
    Q.setFromTriplets(triplets.begin(), triplets.end());
}