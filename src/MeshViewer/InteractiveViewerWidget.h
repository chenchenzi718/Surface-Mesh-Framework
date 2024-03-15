#ifndef INTERACTIVE_VIEWER_WIDGET
#define INTERACTIVE_VIEWER_WIDGET

#include "MeshViewerWidget.h"
#include "../ANN/ANN.h"
#include<queue>
#include<map>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<algorithm>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <set>

/*
	MeshViewerWidget �����˹��� mesh �Ļ�����Ϣ��InteractiveViewerWidget �����Ҫ������仯����
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
	// ���� simplify ��
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
	// ��Ҫʵ���˶����棬�ߣ����ѡȡ�����ƵĹ���
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
	// ���Ȼ������ı�ѡ������ˣ�Ȼ��ִ���˸��� draw_scene
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


#pragma region HomeWork_region
public:
	// ����һ���ṹ��������¼�ߵĺϲ������Լ����ŵ�
	struct EdgeCollapseInfo {
		Mesh::EdgeHandle eh;
		double cost;
		OpenMesh::Vec3d opt_point;
		int count; // ��Щ�߽�ȥ���ȶ���ʱ�Ǹ���ǰ����Ϣ��������
	};
	// ���ñȽϺ���
	struct EdgeCompare {
		bool operator()(const EdgeCollapseInfo& a, const EdgeCollapseInfo& b) {
			// ����MyClass�����е�value���бȽ�
			// Ϊ�˴�����С�ѣ�����Ӧ������true����һ���������ڵڶ�������
			return a.cost > b.cost;
		}
	};


	// ִ�� QEM �㷨
	void QEMSimplifyMesh();

	// �Լ�д�� QEM ���庯��������Ҫ�ǰ� mesh ���޸��� mesh ������Ϊ�����ų����Ͳ��ᱨ��
	void NewMeshSimplification();

	// ���ü򻯵ĵȼ�
	void SetSimplifyLevel(int level)
	{
		// ÿ�� level �仯�������ʱ��Ҫ���ٵ����ٵ㣬Ȼ����в���
		simplify_mesh_vertex =
			(mesh_vertex_num >> level) >= simplify_lowest_vertex_num ? (mesh_vertex_num >> level) : simplify_lowest_vertex_num;
	}

	void Clear(std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare>& edge_que)
	{
		std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare> empty;
		swap(empty, edge_que);
	}

	// ����ÿ�� vertices �� Q ����
	void Cal_Q_ForVertex(Mesh::VertexHandle vertex_h);

	// ����һ���ߺϲ�֮��Ĵ��ۣ��������һ�����ŵ�
	void Cal_Cost_ForEdge(Mesh::EdgeHandle eh);

	// �����µĵ㴦 Q ֵ�Լ��ٽ��ߵ����¼���
	void Update_Local_Variable(Mesh::VertexHandle vertex_h);

	// ֱ�ӽ������� mousePressEvent �ڵĴ���
	bool EdgeCollapse(const EdgeCollapseInfo& edge_collapse_info);


protected:
	int mesh_vertex_num;      // mesh �� vertex ��Ŀ
	int simplify_level = 0;   // �򻯵ȼ�
	int simplify_mesh_vertex; // �򻯺��ж��ٵ�
	double penalty = 1e3;     // ���ñ߽�߳ͷ�
	double cost_min = 1.;  // �����С��ֵ�������⻹��Ͳ�ִ��
	int simplify_lowest_vertex_num = 10;

	std::map<Mesh::VertexHandle, Eigen::Matrix4d> Q2v;
	std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare> edge_que;
	std::map<Mesh::EdgeHandle, int>each_edge_update_num_count;

#pragma endregion
#pragma region other_method_region
public:
	struct AUX_EdgeCollapse
	{
		Mesh::HalfedgeHandle halfeh;
		Mesh::Point newpt;
		Mesh::VertexHandle vto;
		Mesh::VertexHandle vfrom;

		//�жϸ� edge ���ߵ���Ƿ��ѱ����¹��ĵ��ȡ��
		int vto_flag = 0;
		int vfrom_flag = 0;

		Eigen::Matrix4f Q_new;
		float cost;
		bool operator<(const AUX_EdgeCollapse& a) const
		{
			return cost > a.cost;
		}
	};


	void Mesh_Simplification(Mesh& mesh, float ratio);

#pragma endregion
};

#endif