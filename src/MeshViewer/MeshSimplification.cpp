#include "MeshSimplification.h"

MeshSimplification::MeshSimplification(QWidget* parent = 0)
	:InteractiveViewerWidget(parent)
{
	simplify_mesh_vertex = mesh_vertex_num = mesh.n_vertices();

	if (simplify_lowest_vertex_num > mesh_vertex_num) 
		simplify_lowest_vertex_num = 3;
}


MeshSimplification::~MeshSimplification()
{
	// 不需要显示调用父类析构
}

