#ifndef INTERACTIVE_VIEWER_WIDGET
#define INTERACTIVE_VIEWER_WIDGET

#include "MeshViewerWidget.h"
#include "../ANN/ANN.h"

#include"QEMSimplification.h"
#include"SpectralQEM.h"

/*
	MeshViewerWidget 给出了关于 mesh 的基础信息后，InteractiveViewerWidget 完成主要的网格变化工作
*/


class InteractiveViewerWidget : public MeshViewerWidget
{
	Q_OBJECT
public:
	InteractiveViewerWidget(QWidget* parent = 0);
	InteractiveViewerWidget(QGLFormat& _fmt, QWidget* _parent);
	~InteractiveViewerWidget();

	void clearSelectedData()
	{
		selectedVertex.clear();
		selectedFace.clear();
		selectedEdge.clear();
	};

	virtual void clearAllMesh()
	{
		draw_new_mesh = false;
		clearSelectedData();
		MeshViewerWidget::clearAllMesh();
	}

	void setDrawNewMesh(bool draw_new_mesh_)
	{
		draw_new_mesh = draw_new_mesh_;
		updateGL();
	}

	void set_mesh_ref(Mesh& mesh_)
	{
		clearAllMesh();
		mesh = mesh_;
		initMesh(); 
		setDrawMode(FLAT_POINTS);
		setMouseMode(TRANS);
	};

	void edit_undo_viewer()
	{
		--mesh_vector_index;
		mesh = mesh_vector[mesh_vector_index];
		emit set_edit_redo_enable_viewer_signal( true );
		if(mesh_vector_index == 0)
		{
			emit set_edit_undo_enable_viewer_signal( false );
		}
		updateGL();
	}
	void edit_redo_viewer()
	{
		++mesh_vector_index;
		mesh = mesh_vector[mesh_vector_index];
		emit set_edit_undo_enable_viewer_signal( true );
		if(mesh_vector_index == mesh_vector.size() - 1)
		{
			emit set_edit_redo_enable_viewer_signal( false );
		}
		updateGL();
	}

signals:
	void mouse_press_signal(Mesh::Point P);
	void mouse_move_signal(OpenMesh::Vec3d xy);
	void mouse_release_signal(Mesh::Point  P);
	void draw_from_out_signal();

	void setMouseMode_signal(int);

	void set_edit_undo_enable_viewer_signal(bool);
	void set_edit_redo_enable_viewer_signal(bool);

public slots:
	void render_text_slot(OpenMesh::Vec3d pos, QString str);
	void set_t2_mouse_mode(int tm)
	{
		t2_mode_ = tm;
	}

public:
	// 增加 simplify 项
	enum { TRANS, POINTPICK, VERTEXPICK, EDGEPICK, FACEPICK, EDGECOLLAPSE, EDGEFLIP, 
		EDGESPLIT , MOVE, T2_MODE, N_MODE, SIMPLIFY };
	void setMouseMode(int mm);
	int mouseMode() const { return mouse_mode_; }

protected:
	virtual void mousePressEvent(QMouseEvent *_event);
	virtual void mouseReleaseEvent(QMouseEvent *_event);
	virtual void mouseMoveEvent(QMouseEvent *_event);
	virtual void wheelEvent(QWheelEvent* _event);
	int mouse_mode_;
	int t2_mode_;

protected:
	// 主要实现了对于面，线，点的选取并绘制的功能
	void pick_vertex(int x,int y);
	void pick_face(int x,int y);
	void pick_edge(int x,int y);
	void pick_point(int x,int y);
	void move_point_based_lastVertex(int x,int y);

	int find_vertex_using_selected_point();
	int find_face_using_selected_point();
	int find_edge_using_selected_point();

	void buildIndex();
	ANNkd_tree* kdTree;

	void draw_interactive_portion(int drawmode);
	void draw_interactive_portion_mesh2();
	void draw_selected_point();
	void draw_selected_vertex();
	void draw_selected_face();
	void draw_selected_edge();
	// 首先绘出上面的被选择的拓扑，然后执行了父类 draw_scene
	virtual void draw_scene(int drawmode);
	bool draw_new_mesh;

protected:
	double selectedPoint[3];
	std::vector<int> selectedVertex;
	int lastestVertex;
	std::vector<int> selectedFace;
	int lastestFace;
	std::vector<int> selectedEdge;
	int lastestEdge;

protected:
	void dragEnterEvent(QDragEnterEvent *event);
	void dropEvent(QDropEvent *event);

public:

private:

#pragma region Auxiliary_function
public:
	void inverse_mesh_connectivity();
	void scale_mesh_using_BBox(int max_len);//don't change the ratio of the xyz
	void tetgen_save();
	void split_quad_mesh();
	void transform_mesh(const std::vector<double>& m);
	void find_vertex_by_id(int id);
	void find_face_by_id(int id);
	void find_edge_by_id(int id);
	void find_vertex_by_valance(int valance);
	void delete_vertex_valence_four();
	void delete_vertex_valence_three();
	void split_vertex_valence_eight();
#pragma endregion


#pragma region hw1
public:
	QEMSimplification qem;

	void QEMSimplifyMesh();

public slots:
	void QEM_Simplification_Ratio_Changed(int value);

#pragma endregion


#pragma region hw_cubic_stylization
public:
	SpectralQEM sp_qem;

	void SpectralQEMSimplifyMesh();
#pragma endregion
};

#endif