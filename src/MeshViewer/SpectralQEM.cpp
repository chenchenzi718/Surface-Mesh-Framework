#include"SpectralQEM.h"

#include <igl/per_vertex_normals.h>
#include <igl/arap.h>

SpectralQEM::SpectralQEM()
{

}


SpectralQEM::~SpectralQEM()
{

}


void SpectralQEM::ARAP_kernel(Mesh& mesh)
{
    change_openmesh_to_eigen(mesh, this->V, this->F);

    // ���㶥�㷨��
    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);

    // ���ó�ʼ����λ��
    U = V;

    // ������κ�Ŀ��Ƶ�
    b.resize(2, 1); // �����������������Ƶ�
    b << 1, 3;      // Լ����������

    Eigen::MatrixXd bc(2, 3);  // bc ������Ӧ�� b �ĳ�����ͬ
    bc.row(0) = V.row(1) + Eigen::RowVector3d(0.5, 0, 0); // ������ 1 �����ƶ� 0.5 ����λ
    bc.row(1) = V.row(3) + Eigen::RowVector3d(0.3, 0, 0); // ������ 3 �����ƶ� 0.5 ����λ

    // ���� ARAP �㷨����
    igl::ARAPData arap_data;
    arap_data.with_dynamics = false;
    arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;

    // Ԥ���� ARAP
    arap_data.max_iter = 10;
    igl::arap_precomputation(V, F, V.cols(), b, arap_data);

    // ������� ARAP ����
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


