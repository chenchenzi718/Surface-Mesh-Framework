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


#pragma region HomeWork_region
public:
	// 设置一个结构体用来记录边的合并代价以及最优点
	struct EdgeCollapseInfo {
		Mesh::EdgeHandle eh;
		double cost;
		OpenMesh::Vec3d opt_point;
		int count; // 有些边进去优先队列时是更新前的信息，不作数
	};
	// 设置比较函数
	struct EdgeCompare {
		bool operator()(const EdgeCollapseInfo& a, const EdgeCollapseInfo& b) {
			// 根据MyClass对象中的value进行比较
			// 为了创建最小堆，这里应当返回true当第一个参数大于第二个参数
			return a.cost > b.cost;
		}
	};


	// 执行 QEM 算法
	void QEMSimplifyMesh();

	// 自己写的 QEM 主体函数，发现要是把 mesh 的修改与 mesh 另外作为函数放出来就不会报错
	void NewMeshSimplification();

	// 设置简化的等级
	void SetSimplifyLevel(int level)
	{
		// 每次 level 变化，计算此时需要减少到多少点，然后进行操作
		simplify_mesh_vertex =
			(mesh_vertex_num >> level) >= simplify_lowest_vertex_num ? (mesh_vertex_num >> level) : simplify_lowest_vertex_num;
	}

	void Clear(std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare>& edge_que)
	{
		std::priority_queue<EdgeCollapseInfo, std::vector<EdgeCollapseInfo>, EdgeCompare> empty;
		swap(empty, edge_que);
	}

	// 设置每个 vertices 的 Q 矩阵
	void Cal_Q_ForVertex(Mesh::VertexHandle vertex_h);

	// 计算一条边合并之后的代价，并计算出一个最优点
	void Cal_Cost_ForEdge(Mesh::EdgeHandle eh);

	// 更新新的点处 Q 值以及临近边的重新计算
	void Update_Local_Variable(Mesh::VertexHandle vertex_h);

	// 直接借鉴上面的 mousePressEvent 内的处理
	bool EdgeCollapse(const EdgeCollapseInfo& edge_collapse_info);


protected:
	int mesh_vertex_num;      // mesh 的 vertex 数目
	int simplify_level = 0;   // 简化等级
	int simplify_mesh_vertex; // 简化后有多少点
	double penalty = 1e3;     // 设置边界边惩罚
	double cost_min = 1.;  // 如果最小的值弹出比这还大就不执行
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

		//判断该 edge 两边点对是否已被更新过的点对取代
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