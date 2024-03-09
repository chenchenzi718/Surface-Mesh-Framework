#ifndef MESH_SIMPLIFICATION
#define MESH_SIMPLIFICATION


#include "InteractiveViewerWidget.h"


/*
	���� HW1 �ڵ���������̣�ʹ�� QEM �㷨
	����Ҳ���д drawScene()����ô updateGL ����˳����� 
	updateGL()���������ػ津�����ƣ��� QGLViewerWidget::paintGL() �� InteractiveViewerWidget::draw_scene()
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
		// ÿ�� level �仯�������ʱ��Ҫ���ٵ����ٵ㣬Ȼ����в���
		simplify_mesh_vertex =
			(mesh_vertex_num >> level) >= simplify_lowest_vertex_num ? (mesh_vertex_num >> level) : simplify_lowest_vertex_num;
	}




private:
	// mesh �� vertex ��Ŀ
	int mesh_vertex_num;

	// �򻯵ȼ�
	int simplify_level = 0;
	// �򻯺��ж��ٵ�
	int simplify_mesh_vertex;

	int simplify_lowest_vertex_num = 10;
};

#endif