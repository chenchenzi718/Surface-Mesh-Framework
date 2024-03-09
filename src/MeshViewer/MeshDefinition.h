#ifndef MESHDEFINITION_H
#define MESHDEFINITION_H

#include <OpenMesh/Core/Geometry/VectorT.hh>
//#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

/*
	制定网格的属性
	定义网格特征（MeshTraits）：
		通过继承OpenMesh的DefaultTraits，自定义了一个MeshTraits结构体，用于定制网格的属性。这里包括：
			点（Point）和法线（Normal）使用OpenMesh::Vec3d类型，即三维向量。
			启用顶点、面、边和半边的状态属性，这些状态可用于标记元素是否被删除、锁定等。
			定义了额外的特征（Traits）和私有成员变量用于半边、边和顶点，如face_he_var（可能用于存储与面相关的半边信息）和new_pos（用于存储顶点的新位置，可能用于网格变形或参数化）等，以及相关的公开接口函数。
			定义了网格类型：使用MeshTraits作为模板参数，通过typedef定义了Mesh类型为OpenMesh::PolyMesh_ArrayKernelT<MeshTraits>。这表明这个头文件中将处理的主要数据结构是多边形网格。
*/

struct MeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;

	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);

	FaceTraits
	{
	};

	EdgeTraits
	{
	};

	HalfedgeTraits
	{
		HalfedgeT() :face_he_var(-1)
		{
		};
	private:
		int face_he_var;
	public:
		int get_face_he_var() { return face_he_var; };
		void set_face_he_var(int fhe) { face_he_var = fhe; };
	};

	VertexTraits
	{
		VertexT() : new_pos_fixed(false)
		{
		};
	private:
		OpenMesh::Vec3d new_pos;//can be used for deformation and parameterization
		bool new_pos_fixed;
	public:
		void set_New_Pos(const OpenMesh::Vec3d& n_p){ new_pos = n_p; };
		OpenMesh::Vec3d& get_New_Pos(){ return new_pos; };
		void set_new_pos_fixed(bool f){ new_pos_fixed = f; };
		bool get_new_pos_fixed(){ return new_pos_fixed; };
	};

};

//typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> Mesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> Mesh;

bool is_flip_ok_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);//just copy the code from openmesh
bool flip_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);

bool check_in_triangle_face(const std::vector<OpenMesh::Vec3d>& tri, const OpenMesh::Vec3d& p);


#endif