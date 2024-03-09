#ifndef MESH_SIMPLIFICATION
#define MESH_SIMPLIFICATION


#include "InteractiveViewerWidget.h"


/*
	处理 HW1 内的网格简化流程，使用 QEM 算法
	如果我不重写 drawScene()，那么 updateGL 调用顺序就是 
	updateGL()（或其他重绘触发机制）→ QGLViewerWidget::paintGL() → InteractiveViewerWidget::draw_scene()
*/
class MeshSimplification : public InteractiveViewerWidget
{
	Q_OBJECT

public:
	MeshSimplification(QWidget* parent = 0);
	MeshSimplification(QGLFormat& _fmt, QWidget* _parent);
	~MeshSimplification();

	void SetSimplifyLevel(int level)
	{
		// 每次 level 变化，计算此时需要减少到多少点，然后进行操作
		simplify_mesh_vertex =
			(mesh_vertex_num >> level) >= simplify_lowest_vertex_num ? (mesh_vertex_num >> level) : simplify_lowest_vertex_num;
	}




private:
	// mesh 的 vertex 数目
	int mesh_vertex_num;

	// 简化等级
	int simplify_level = 0;
	// 简化后有多少点
	int simplify_mesh_vertex;

	int simplify_lowest_vertex_num = 10;
};

#endif