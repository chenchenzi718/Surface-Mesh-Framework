#include"QEMSimplification.h"

QEMSimplification::QEMSimplification(float ratio)
	:ratio(ratio)
{
	penalty = 1e3;
	cost_min = 1.;
	simplify_lowest_vertex_num = 3;
}


QEMSimplification::~QEMSimplification()
{
	Q2v.clear();
	Clear(edge_que); // ������� handle ����ɾ���garbage �ռ�֮��ͱ���
	each_edge_update_num_count.clear();
}


void QEMSimplification::Cal_Q_ForVertex(const Mesh& mesh, Mesh::VertexHandle vertex_h)
{
	// ��ȡ����λ��
	OpenMesh::Vec3d pos = mesh.point(vertex_h);


	Eigen::Matrix4d Q_temp;  // ������¼ÿ��������� Q ����������
	Q_temp.setZero();

	Eigen::Vector4d p;   // ������¼�����ڵ� plane
	p.setZero();


	/*
	// �������㲻�Ǳ߽��
	if (!mesh.is_boundary(vertex_h))
	{
		// ��������
		for (Mesh::ConstVertexFaceIter vf_it = mesh.cvf_iter(vertex_h); vf_it.is_valid(); ++vf_it)
		{
			Mesh::Normal face_normal = mesh.normal(vf_it);
			double d = -face_normal | pos;  // ������
			p = { face_normal[0],face_normal[1],face_normal[2],d };
			Q_temp += p * p.transpose();
		}
		Q2v[vertex_h] = Q_temp;
	}
	else
	{
		// ��ʱΪ�߽�㣬�ù۲���Ƿ�Ϊ�߽�ߣ����������趨�ı߽�߻������������Ǳ߽�㣬������
		// ���������������������ǴӰ�߳�������Ѱ��
		for (Mesh::ConstVertexOHalfedgeIter voh_it = mesh.cvoh_iter(vertex_h); voh_it.is_valid(); ++voh_it)
		{
			Mesh::HalfedgeHandle heh = *voh_it;       // ȡ�����
			Mesh::EdgeHandle eh = mesh.edge_handle(heh);
			Mesh::VertexHandle other_vh = mesh.to_vertex_handle(heh); // ȡ����һ����

			if (!mesh.is_boundary(heh))
			{
				// �����߲��Ǳ߽磬����һ�������Ӧ������� Q_temp
				Mesh::FaceHandle fh = mesh.face_handle(heh);
				Mesh::Normal face_normal = mesh.normal(fh);
				double d = -face_normal | pos;  // ������
				p = { face_normal[0],face_normal[1],face_normal[2],d };
				Q_temp += p * p.transpose();
			}
			else
			{
				// �������Ǳ߽�,��������ߴ�����һ���޴��Լ���棬һ���㴦���������µ�Լ����
				Mesh::HalfedgeHandle prev_heh = mesh.prev_halfedge_handle(heh);

				// ��ð�ߵĶԱ�������ķ���
				Mesh::FaceHandle fheh = mesh.face_handle(mesh.opposite_halfedge_handle(heh));
				Mesh::FaceHandle f_prevheh = mesh.face_handle(mesh.opposite_halfedge_handle(prev_heh));
				Mesh::Normal fheh_normal = mesh.normal(fheh);
				Mesh::Normal f_prevheh_normal = mesh.normal(f_prevheh);

				// ���������߽��ߵķ�������
				OpenMesh::Vec3d heh_tang = mesh.point(other_vh) - mesh.point(vertex_h);
				Mesh::VertexHandle prev_heh_from_v = mesh.from_vertex_handle(prev_heh);
				Mesh::VertexHandle prev_heh_to_v = mesh.to_vertex_handle(prev_heh);
				OpenMesh::Vec3d prev_heh_tang = mesh.point(prev_heh_to_v) - mesh.point(prev_heh_from_v);

				// ���������������ķ���
				OpenMesh::Vec3d virtual_face_normal = OpenMesh::cross(heh_tang, fheh_normal).normalized();
				OpenMesh::Vec3d virtual_prev_face_normal = OpenMesh::cross(prev_heh_tang, f_prevheh_normal).normalized();

				double d = -virtual_face_normal | pos;
				p = { virtual_face_normal[0],virtual_face_normal[1],virtual_face_normal[2],d };
				Q_temp += (p * p.transpose());

				d = -virtual_prev_face_normal | pos;
				p = { virtual_prev_face_normal[0],virtual_prev_face_normal[1],virtual_prev_face_normal[2],d };
				Q_temp += (p * p.transpose());
			}

			// �����������������б߲��Ǳ߽磬�����������Ǳ߽�����
			if (!mesh.is_boundary(eh) && mesh.is_boundary(other_vh))
			{
				// Q_temp *= penalty;
			}
		}
		Q2v[vertex_h] = Q_temp;
	}
	*/

	for (Mesh::ConstVertexFaceIter vf_it = mesh.cvf_iter(vertex_h); vf_it.is_valid(); ++vf_it)
	{
		Mesh::Normal face_normal = mesh.normal(*vf_it);
		double d = -face_normal | pos;  // ������
		p = { face_normal[0],face_normal[1],face_normal[2],d };
		Q_temp += p * p.transpose();
	}
	Q2v[vertex_h] = Q_temp;
}


void QEMSimplification::Cal_Cost_ForEdge(const Mesh& mesh, Mesh::EdgeHandle eh)
{
	EdgeCollapseInfo info;
	Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh);
	Mesh::VertexHandle from_v = mesh.from_vertex_handle(heh);
	Mesh::VertexHandle to_v = mesh.to_vertex_handle(heh);


	Eigen::Matrix4d Q_edge = Q2v[from_v] + Q2v[to_v];
	Eigen::Matrix4d A = Q_edge;
	A(3, 0) = 0.0, A(3, 1) = 0.0, A(3, 2) = 0.0, A(3, 3) = 1.0;

	Eigen::Vector4d b = { 0.0,0.0,0.0,1.0 };

	Eigen::Vector4d cal_vec = A.colPivHouseholderQr().solve(b);

	OpenMesh::Vec3d opt_point(cal_vec[0], cal_vec[1], cal_vec[2]);

	info.eh = eh;
	info.cost = cal_vec.transpose() * Q_edge * cal_vec;
	info.opt_point = opt_point;
	info.count = each_edge_update_num_count[eh];
	edge_que.push(info);
}


void QEMSimplification::Update_Local_Variable(Mesh& mesh, Mesh::VertexHandle vh)
{
	for (auto ve_it = mesh.ve_iter(vh); ve_it.is_valid(); ++ve_it)
	{
		each_edge_update_num_count[*ve_it]++;
		Cal_Cost_ForEdge(mesh, *ve_it);
	}
}


bool QEMSimplification::EdgeCollapse(Mesh& mesh, const EdgeCollapseInfo& edge_collapse_info)
{
	auto eh = edge_collapse_info.eh;
	Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh);
	Mesh::HalfedgeHandle heh_other = mesh.opposite_halfedge_handle(heh);

	bool collapse_ok = false;

	auto from_p = mesh.from_vertex_handle(heh);
	auto to_p = mesh.to_vertex_handle(heh);

	auto sp = edge_collapse_info.opt_point;
	Mesh::VertexHandle v_result;

	auto Q_result = Q2v[from_p] + Q2v[to_p];


	if (mesh.is_collapse_ok(heh))
	{
		mesh.set_point(to_p, sp);
		mesh.collapse(heh);
		collapse_ok = true;
		v_result = to_p;
	}
	else if (mesh.is_collapse_ok(heh_other))
	{
		mesh.set_point(from_p, sp);
		mesh.collapse(heh_other);
		collapse_ok = true;
		v_result = from_p;
	}
	else
	{
		return collapse_ok;
	}


	Q2v[v_result] = Q_result;
	Update_Local_Variable(mesh, v_result);


	return collapse_ok;
}


void QEMSimplification::QEMSimplifyMesh(Mesh& mesh)
{
	uint mesh_vertex_num = mesh.n_vertices();
	uint simplify_mesh_vertex = static_cast<uint>(ratio * mesh_vertex_num);
	if (simplify_mesh_vertex < simplify_lowest_vertex_num)
		simplify_mesh_vertex = simplify_lowest_vertex_num;

	// ���ȼ��� mesh ���е㴦�� Q ֵ���� Q2v ��
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		Cal_Q_ForVertex(mesh, *v_it);
	}

	// ����ÿ���ߵ� cost
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		each_edge_update_num_count[*e_it] = 0;
		Cal_Cost_ForEdge(mesh, *e_it);
	}

	int collapse_num = 0;
	// ����ѭ��
	while (collapse_num < (mesh_vertex_num - simplify_mesh_vertex) && !edge_que.empty())
	{
		// ȡ����С�� cost ��
		EdgeCollapseInfo edgeinfo = edge_que.top();
		edge_que.pop();

		auto eh = edgeinfo.eh;
		if (mesh.status(eh).deleted())
			continue;
		if (edgeinfo.count != each_edge_update_num_count[eh])
			continue;

		EdgeCollapse(mesh, edgeinfo);
		collapse_num++;
	}

	Q2v.clear();
	Clear(edge_que); // ������� handle ����ɾ���garbage �ռ�֮��ͱ���
	each_edge_update_num_count.clear();

	mesh.garbage_collection();
}