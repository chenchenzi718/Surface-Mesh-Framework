#include <QMouseEvent>
#include <QLineEdit>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QtCore>
#include <QUrl>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "InteractiveViewerWidget.h"

InteractiveViewerWidget::InteractiveViewerWidget(QWidget* parent /* = 0 */)
	:MeshViewerWidget(parent)
{
	draw_new_mesh = false;
	clearSelectedData();
	kdTree = NULL;

	simplify_mesh_vertex = mesh_vertex_num = mesh.n_vertices();

	if (simplify_lowest_vertex_num > mesh_vertex_num)
		simplify_lowest_vertex_num = 3;
}

InteractiveViewerWidget::InteractiveViewerWidget(QGLFormat& _fmt, QWidget* _parent)
:MeshViewerWidget(_fmt, _parent)
{
	draw_new_mesh = false;
	clearSelectedData();
	kdTree = NULL;

	simplify_mesh_vertex = mesh_vertex_num = mesh.n_vertices();

	if (simplify_lowest_vertex_num > mesh_vertex_num)
		simplify_lowest_vertex_num = 3;
}

InteractiveViewerWidget::~InteractiveViewerWidget()
{
	if(kdTree) delete kdTree;
}

void InteractiveViewerWidget::setMouseMode(int mm)
{
	if(mouse_mode_ != T2_MODE)
	{
		mouse_mode_ = mm;
		if( TRANS != mouse_mode_ )
		{ buildIndex(); }

		// 如果 mousemode 变成了 SIMPLIFY，那么就进行一波网格简化操作
		if(SIMPLIFY == mouse_mode_)
		{
			this->QEMSimplifyMesh();
		}

		emit setMouseMode_signal(mm);
	}
}

void InteractiveViewerWidget::mousePressEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mousePressEvent(_event);
	}
	else
	{
		if(mouse_mode_ != T2_MODE)
		{
			pick_point( _event->x(), _event->y() );
			if(mouse_mode_ == VERTEXPICK)
			{
				pick_vertex( _event->x(), _event->y() );
			}
			else if(mouse_mode_ == FACEPICK)
			{
				pick_face( _event->x(), _event->y() );
			}
			else if(mouse_mode_ == EDGEPICK)
			{
				pick_edge( _event->x(), _event->y() );
			}
			else if(mouse_mode_ == POINTPICK)
			{
			}
			else if( mouse_mode_ == MOVE )
			{
				pick_vertex( _event->x(), _event->y() );//set the selected handle
			}
			else if(mouse_mode_ == EDGECOLLAPSE)
			{
				int desired_edge = find_edge_using_selected_point();
				if(desired_edge >= 0) 
				{
					Mesh::HalfedgeHandle heh = mesh.halfedge_handle( mesh.edge_handle(desired_edge), 0 );
					OpenMesh::Vec3d from_p = mesh.point(mesh.from_vertex_handle(heh));
					OpenMesh::Vec3d to_p = mesh.point(mesh.to_vertex_handle(heh));
					OpenMesh::Vec3d sp(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
					bool collapse_ok = true;
					if( (sp-from_p).sqrnorm() > (to_p-sp).sqrnorm() )
					{
						if( mesh.is_collapse_ok(heh) )
						{
							mesh.collapse(heh);
						}
						else
						{
							collapse_ok = false;
							printf("[%d] Collapse Not OK!\n", desired_edge);
						}
					}
					else
					{
						heh = mesh.opposite_halfedge_handle(heh);
						if( mesh.is_collapse_ok(heh) )
						{
							mesh.collapse(heh);
						}
						else
						{
							collapse_ok = false;
							printf("[%d] Collapse Not OK!\n", desired_edge);
						}
					}
					if(collapse_ok)
					{
						mesh.garbage_collection();
						buildIndex();
						if( mesh_vector.size() - 1 > mesh_vector_index )
						{
							mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
						}
						mesh_vector.push_back( mesh ); mesh_vector_index += 1;
						emit set_edit_undo_enable_viewer_signal( true );
						emit set_edit_redo_enable_viewer_signal( false );
					}
					clearSelectedData();
				}
			}
			else if (mouse_mode_ == EDGEFLIP)
			{
				int desired_edge = find_edge_using_selected_point();
				if(desired_edge >= 0) 
				{
					Mesh::EdgeHandle eh = mesh.edge_handle(desired_edge);
					if( is_flip_ok_openmesh(eh, mesh))
					{
						flip_openmesh(eh, mesh);
						if( mesh_vector.size() - 1 > mesh_vector_index )
						{
							mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
						}
						mesh_vector.push_back( mesh ); mesh_vector_index += 1;
						emit set_edit_undo_enable_viewer_signal( true );
						emit set_edit_redo_enable_viewer_signal( false );
					}
					else
					{
						printf("[%d] Flip Not OK!\n", desired_edge);
					}
					clearSelectedData();
				}
			}
			else if (mouse_mode_ == EDGESPLIT)
			{
				int desired_edge = find_edge_using_selected_point();
				if(desired_edge >= 0) 
				{
					Mesh::EdgeHandle eh = mesh.edge_handle(desired_edge);
					Mesh::HalfedgeHandle heh = mesh.halfedge_handle( eh, 0 );
					Mesh::HalfedgeHandle heh_ = mesh.halfedge_handle( eh, 1 );
					Mesh::VertexHandle vh0 = mesh.to_vertex_handle(heh);
					Mesh::VertexHandle vh1 = mesh.to_vertex_handle(heh_);
					OpenMesh::Vec3d s = mesh.point( vh1 );
					OpenMesh::Vec3d e = mesh.point( vh0 );
					Mesh::VertexHandle vh = mesh.add_vertex( (s + e)*0.5 );
					std::vector<Mesh::VertexHandle> one_face(3);
					if(mesh.is_boundary(eh))
					{
						if(Mesh::InvalidFaceHandle != mesh.face_handle(heh))
						{
							Mesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
							mesh.delete_edge(eh, false); mesh.garbage_collection();
							one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh.add_face(one_face);
							one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh.add_face(one_face);
						}
						else
						{
							Mesh::VertexHandle vh3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh_));
							mesh.delete_edge(eh, false); mesh.garbage_collection();
							one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh.add_face(one_face);
							one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh.add_face(one_face);
						}
					}
					else
					{
						Mesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
						Mesh::VertexHandle vh3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh_));
						mesh.delete_edge(eh, false); mesh.garbage_collection();
						one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh.add_face(one_face);
						one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh.add_face(one_face);
						one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh.add_face(one_face);
						one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh.add_face(one_face);
					}

					mesh.update_normals();
					buildIndex();
					clearSelectedData();

					if( mesh_vector.size() - 1 > mesh_vector_index )
					{
						mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
					}
					mesh_vector.push_back( mesh ); mesh_vector_index += 1;
					emit set_edit_undo_enable_viewer_signal( true );
					emit set_edit_redo_enable_viewer_signal( false );
				}
			}
		}
	}
	updateGL();
}

void InteractiveViewerWidget::mouseMoveEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mouseMoveEvent(_event);
	}
	else
	{
		if( mouse_mode_ != T2_MODE)
		{
			if( mouse_mode_ == MOVE )
			{
				move_point_based_lastVertex( _event->x(), _event->y() );
				Mesh::Point P(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
				mesh.set_point( mesh.vertex_handle(lastestVertex), P );
				updateGL();
			}
		}
		else
		{
			
		}
		
	}
}

void InteractiveViewerWidget::mouseReleaseEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mouseMoveEvent(_event);
	}
	else
	{
		if(mouse_mode_ != T2_MODE )
		{
			if( mouse_mode_ == MOVE )
			{
				move_point_based_lastVertex( _event->x(), _event->y() );
				Mesh::Point P(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
				mesh.set_point( mesh.vertex_handle(lastestVertex), P );
				selectedVertex.clear();

				// 修改原框架此处代码
				mesh_vector.push_back(mesh); mesh_vector_index += 1;
				emit set_edit_undo_enable_viewer_signal(true);
				emit set_edit_redo_enable_viewer_signal(false);

				updateGL();
			}
		}
		else
		{
		}
	}
	
}

void InteractiveViewerWidget::wheelEvent(QWheelEvent* _event)
{
	if(mouse_mode_ != N_MODE && mouse_mode_ != T2_MODE)
	{
		MeshViewerWidget::wheelEvent(_event);
	}
}

void InteractiveViewerWidget::dragEnterEvent(QDragEnterEvent* event)
{
	if( event->mimeData()->hasFormat("text/uri-list") )
	{
		event->acceptProposedAction();
	}
}

void InteractiveViewerWidget::dropEvent(QDropEvent* event)
{
	QList<QUrl> urls = event->mimeData()->urls();
	if( urls.isEmpty() )
		return;
	QString fileName = urls.first().toLocalFile();
	if (fileName.isEmpty())
		return;

	if( fileName.endsWith(".off") || fileName.endsWith(".obj") || fileName.endsWith(".stl") || fileName.endsWith(".ply"))
	{
		if( openMesh(fileName.toLocal8Bit()))
		{
			emit(loadMeshOK(true,fileName));
			setDrawMode(FLAT_POINTS);
			setMouseMode(TRANS);
		}
		else
		{
			emit(loadMeshOK(false,"No Mesh"));
		}
	}
}

void InteractiveViewerWidget::pick_vertex(int x,int y)
{
	int r = find_vertex_using_selected_point();
	lastestVertex = r;
	printf("Select Vertex : %d\n", r);
	std::vector<int>::iterator it;
	if( (it = std::find(selectedVertex.begin(),selectedVertex.end(), r)) == selectedVertex.end() )
	{
		selectedVertex.push_back(r);
	}
	else
	{
		selectedVertex.erase(it);
	}

	updateGL();
}
void InteractiveViewerWidget::pick_face(int x,int y)
{
	int desiredFace = find_face_using_selected_point();
	if(desiredFace < 0) return;
	lastestFace = desiredFace;
	printf("Select Face : %d\n", desiredFace);
	std::vector<int>::iterator it;
	if( (it = std::find(selectedFace.begin(),selectedFace.end(),desiredFace)) == selectedFace.end() )
	{
		selectedFace.push_back(desiredFace);
	}
	else
	{
		selectedFace.erase(it);
	}
	updateGL();
}
void InteractiveViewerWidget::pick_edge(int x,int y)
{
	int desiredEdge = find_edge_using_selected_point();
	if(desiredEdge < 0) return;
	lastestEdge = desiredEdge;
	printf("Select Edge : %d\n", desiredEdge);
	std::vector<int>::iterator it;
	if( (it = std::find(selectedEdge.begin(),selectedEdge.end(),desiredEdge)) == selectedEdge.end() )
	{
		selectedEdge.push_back(desiredEdge);
	}
	else
	{
		selectedEdge.erase(it);
	}
	updateGL();
}
void InteractiveViewerWidget::pick_point(int x,int y)
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = double(x);
	GLdouble winY = double( height() - y );
	GLfloat winZ = 0.0;
	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
	gluUnProject(winX, winY, (GLdouble)winZ, &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}

void InteractiveViewerWidget::move_point_based_lastVertex(int x,int y)
{
	if(lastestVertex<0 || lastestVertex>=mesh.n_vertices())
	{
		return;
	}
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = 0.0;
	GLdouble winY = 0.0;
	GLdouble winZ = 0.0;
	OpenMesh::Vec3d p = mesh.point(mesh.vertex_handle(lastestVertex));
	gluProject(p[0], p[1], p[2],  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &winX, &winY, &winZ);
	
	gluUnProject((GLdouble)(x), (GLdouble)( height() - y ), winZ,  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}

int InteractiveViewerWidget::find_vertex_using_selected_point()
{
	ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	kdTree->annkSearch(tp, 1, nnIdx, dists);
	return nnIdx[0];
}

int InteractiveViewerWidget::find_face_using_selected_point()
{
	int rv = find_vertex_using_selected_point();
	Mesh::VertexFaceIter vf_it = mesh.vf_iter( mesh.vertex_handle(rv) );
	int desiredFace = -1; //double minLen = 10*radius();
	// 输入的不是三角网格时出错，Vec3d 超出数组
	std::vector<OpenMesh::Vec3d> tri_p(3); int tri_count = 0;
	Mesh::Point resultP(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
	for( vf_it; vf_it; ++vf_it )
	{
		tri_count = 0;
		for(Mesh::FaceVertexIter fv_it = mesh.fv_iter(vf_it.handle()); fv_it; ++fv_it)
		{
			tri_p[tri_count] = mesh.point(fv_it); ++tri_count;
		}
		if( check_in_triangle_face(tri_p, resultP) )
		{
			desiredFace = vf_it.handle().idx(); break;
		}
	}
	if(desiredFace < 0)
	{
		for(Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
		{
			tri_count = 0;
			for(Mesh::FaceVertexIter fv_it = mesh.fv_iter(f_it.handle()); fv_it; ++fv_it)
			{
				tri_p[tri_count] = mesh.point(fv_it); ++tri_count;
			}
			if( check_in_triangle_face(tri_p, resultP) )
			{
				desiredFace = f_it.handle().idx(); break;
			}
		}
	}

	return  desiredFace;
}

int InteractiveViewerWidget::find_edge_using_selected_point()
{
	int desiredFace = find_face_using_selected_point(); if(desiredFace < 0) return -1;
	Mesh::FaceHandle fh = mesh.face_handle(desiredFace);
	double min_len= 1e30; int desiredEdge = -1;
	Mesh::Point resultP(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
	for(Mesh::FaceHalfedgeIter fhe_it = mesh.fh_iter(fh); fhe_it; ++fhe_it)
	{
		OpenMesh::Vec3d s = mesh.point( mesh.from_vertex_handle(fhe_it) );
		OpenMesh::Vec3d e = mesh.point( mesh.to_vertex_handle(fhe_it) );
		double dis = OpenMesh::cross(resultP - s, resultP - e).norm() / (s - e).norm();
		if(dis < min_len){ min_len = dis; desiredEdge = mesh.edge_handle(fhe_it.handle()).idx(); }
	}
	
	return desiredEdge;
}

// 这段代码的目的是根据一个3D网格的顶点数据构建一个k-d树，以便于高效地执行最近邻查询
void InteractiveViewerWidget::buildIndex()
{
	if(mesh.n_vertices() == 0)
		return;

	Mesh::VertexIter v_it(mesh.vertices_begin());
	Mesh::VertexIter v_end(mesh.vertices_end());
	Mesh::Point p;
	unsigned nv = mesh.n_vertices();
	ANNpointArray dataPts = annAllocPts(nv, 3);
	int count = 0;
	for(; v_it != v_end; ++v_it)
	{
		p = mesh.point(v_it);
		dataPts[count][0] = p[0]; dataPts[count][1] = p[1]; dataPts[count][2] = p[2];
		++count;
	}

	if(kdTree) delete kdTree;
	kdTree = new ANNkd_tree(dataPts, nv, 3);
}

//with the first mesh
void InteractiveViewerWidget::draw_interactive_portion(int drawmode)
{
	glViewport ( 0,0, width(),height());
	glMatrixMode( GL_PROJECTION );
	glLoadMatrixd( &ProjectionMatrix[0] );
	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( &ModelViewMatrix[0] );

	
	emit draw_from_out_signal();

	{
		//draw select vertex, face, edge.
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

		glPointSize(1);

		switch(mouse_mode_)
		{
		case POINTPICK:
			draw_selected_point();
			break;
		case VERTEXPICK:
			draw_selected_vertex();
			break;
		case FACEPICK:
			draw_selected_face();
			break;
		case EDGEPICK:
			draw_selected_edge();
			break;
		default:
			draw_selected_vertex();
			draw_selected_face();
			draw_selected_edge();
			break;
		}
	}

	if(draw_new_mesh)
	{
		draw_scene_mesh(drawmode);
	}
}

//with the second mesh
void InteractiveViewerWidget::draw_interactive_portion_mesh2()
{
	return;
}

void InteractiveViewerWidget::draw_selected_point()
{
	glColor3f(1.0, 0.5, 0.0);
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3d(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
	glEnd();
	glPointSize(1);
}

void InteractiveViewerWidget::draw_selected_vertex()
{
	if( selectedVertex.size() > 0 )
	{
		Mesh::Point p;
		glColor3f(1.0, 0.5, 0.0);
		glPointSize(12);
		glBegin(GL_POINTS);
		for(unsigned int i=0;i<selectedVertex.size();++i)
		{
			p = mesh.point( mesh.vertex_handle(selectedVertex[i]) );
			glVertex3dv(p.data());
		}
		glEnd();
		glPointSize(1);
	}
}

void InteractiveViewerWidget::draw_selected_face()
{
	if( selectedFace.size() > 0 )
	{
		glColor3f(1.0, 0.5, 1.0);
		Mesh::Point p;
		Mesh::ConstFaceVertexIter fv_it;
		Mesh::FaceHandle f_handle;
		for( unsigned int i=0; i<selectedFace.size(); ++i )
		{
			f_handle = mesh.face_handle(selectedFace[i]);
			fv_it = mesh.fv_iter(f_handle);
			glBegin(GL_POLYGON);
			for( fv_it; fv_it; ++fv_it )
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
			glEnd();
		}
	}
}

void InteractiveViewerWidget::draw_selected_edge()
{
	if( selectedEdge.size() > 0)
	{
		glColor3f(1.0, 0.5, 1.0);
		Mesh::Point p1; Mesh::Point p2;
		Mesh::EdgeHandle e_handle;
		Mesh::HalfedgeHandle he_handle;
		for(unsigned int i=0;i<selectedEdge.size();++i)
		{
			e_handle = mesh.edge_handle(selectedEdge[i]);
			he_handle = mesh.halfedge_handle( e_handle, 0 );
			p1 = mesh.point( mesh.from_vertex_handle( he_handle ) );
			p2 = mesh.point( mesh.to_vertex_handle( he_handle ) );
			glBegin(GL_LINES);
			glVertex3dv( p1.data() );
			glVertex3dv( p2.data() );
			glEnd();
		}
	}
}

void InteractiveViewerWidget::draw_scene(int drawmode)
{
	if (!mesh.n_vertices()) { return; }
	draw_interactive_portion_mesh2();
	draw_interactive_portion(drawmode);

	if( !draw_new_mesh )
	{
		MeshViewerWidget::draw_scene(drawmode);
	}
}

void InteractiveViewerWidget::render_text_slot(OpenMesh::Vec3d pos, QString str)
{
	/*GLdouble  winX, winY, winZ;
	GLint     viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	gluProject(pos[0],pos[1],pos[2],&ModelViewMatrix[0][0],&ProjectionMatrix[0][0],viewport,&winX,&winY,&winZ);
	int x = (long)winX;
	int y = viewport[3]-(long)winY;
	render_text(x,y,str);*/
	render_text(pos[0],pos[1],pos[2],str);
}

void InteractiveViewerWidget::Cal_Q_ForVertex(Mesh::VertexHandle vertex_h)
{
	// 获取顶点位置
	OpenMesh::Vec3d pos = mesh.point(vertex_h);


	Eigen::Matrix4d Q_temp;  // 用来记录每个面带来的 Q 矩阵，最后叠加
	Q_temp.setZero();

	Eigen::Vector4d p;   // 用来记录文章内的 plane
	p.setZero();

	/*
	// 如果这个点不是边界点
	if (!mesh.is_boundary(vertex_h))
	{
		// 正常计算
		for (Mesh::ConstVertexFaceIter vf_it = mesh.cvf_iter(vertex_h); vf_it.is_valid(); ++vf_it)
		{
			Mesh::Normal face_normal = mesh.normal(vf_it);
			double d = -face_normal | pos;  // 计算点乘
			p = { face_normal[0],face_normal[1],face_normal[2],d };
			Q_temp += p * p.transpose();
		}
		Q2v[vertex_h] = Q_temp;
	}
	else
	{
		// 此时为边界点，得观察边是否为边界边，这里我们设定的边界边还包括两个点是边界点，但边有
		// 两个邻面的情况，这里我们从半边出发进行寻找
		for (Mesh::ConstVertexOHalfedgeIter voh_it = mesh.cvoh_iter(vertex_h); voh_it.is_valid(); ++voh_it)
		{
			Mesh::HalfedgeHandle heh = *voh_it;       // 取出半边 
			Mesh::EdgeHandle eh = mesh.edge_handle(heh);
			Mesh::VertexHandle other_vh = mesh.to_vertex_handle(heh); // 取出另一个点

			if (!mesh.is_boundary(heh))
			{
				// 如果半边不是边界，计算一次这个对应面带来的 Q_temp
				Mesh::FaceHandle fh = mesh.face_handle(heh);
				Mesh::Normal face_normal = mesh.normal(fh);
				double d = -face_normal | pos;  // 计算点乘
				p = { face_normal[0],face_normal[1],face_normal[2],d };
				Q_temp += p * p.transpose();
			}
			else
			{
				// 如果半边是边界,在这条半边处加上一个巨大的约束面，一个点处会有两个新的约束面
				Mesh::HalfedgeHandle prev_heh = mesh.prev_halfedge_handle(heh);
				
				// 获得半边的对边所在面的法向
				Mesh::FaceHandle fheh = mesh.face_handle(mesh.opposite_halfedge_handle(heh)); 
				Mesh::FaceHandle f_prevheh = mesh.face_handle(mesh.opposite_halfedge_handle(prev_heh));
				Mesh::Normal fheh_normal = mesh.normal(fheh);
				Mesh::Normal f_prevheh_normal = mesh.normal(f_prevheh);

				// 计算两个边界半边的方向向量
				OpenMesh::Vec3d heh_tang = mesh.point(other_vh) - mesh.point(vertex_h);
				Mesh::VertexHandle prev_heh_from_v = mesh.from_vertex_handle(prev_heh);
				Mesh::VertexHandle prev_heh_to_v = mesh.to_vertex_handle(prev_heh);
				OpenMesh::Vec3d prev_heh_tang = mesh.point(prev_heh_to_v) - mesh.point(prev_heh_from_v);

				// 获得两个虚拟邻面的法向
				OpenMesh::Vec3d virtual_face_normal = OpenMesh::cross(heh_tang, fheh_normal).normalized();
				OpenMesh::Vec3d virtual_prev_face_normal = OpenMesh::cross(prev_heh_tang, f_prevheh_normal).normalized();
				
				double d = -virtual_face_normal | pos;
				p = { virtual_face_normal[0],virtual_face_normal[1],virtual_face_normal[2],d };
				Q_temp += (p * p.transpose());

				d = -virtual_prev_face_normal | pos;
				p = { virtual_prev_face_normal[0],virtual_prev_face_normal[1],virtual_prev_face_normal[2],d };
				Q_temp += (p * p.transpose());
			}

			// 除了上面的情况，还有边不是边界，但是两个点是边界的情况
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
		Mesh::Normal face_normal = mesh.normal(vf_it);
		double d = -face_normal | pos;  // 计算点乘
		p = { face_normal[0],face_normal[1],face_normal[2],d };
		Q_temp += p * p.transpose();
	}
	Q2v[vertex_h] = Q_temp;
}


void InteractiveViewerWidget::Cal_Cost_ForEdge(Mesh::EdgeHandle eh)
{
	EdgeCollapseInfo info;
	Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh);
	Mesh::VertexHandle from_v = mesh.from_vertex_handle(heh);
	Mesh::VertexHandle to_v = mesh.to_vertex_handle(heh);


	Eigen::Matrix4d Q_edge = Q2v[from_v] + Q2v[to_v];
	Eigen::Matrix4d Q_cal = Q_edge;
	Q_cal(3, 0) = 0.0, Q_cal(3, 1) = 0.0, Q_cal(3, 2) = 0.0, Q_cal(3, 3) = 1.0;

	OpenMesh::Vec3d opt_point;
	Eigen::Vector4d cal_vec;
	Eigen::Vector4d cal_b = { 0.0,0.0,0.0,1.0 };

	cal_vec = Q_cal.colPivHouseholderQr().solve(cal_b);

	opt_point = { cal_vec[0], cal_vec[1], cal_vec[2] };

	info.eh = eh;
	info.cost = cal_vec.transpose() * Q_edge * cal_vec;
	info.opt_point = opt_point;
	info.count = each_edge_update_num_count[eh];
	edge_que.push(info);
}


void InteractiveViewerWidget::Update_Local_Variable(Mesh::VertexHandle vh)
{
	for (auto ve_it = mesh.ve_iter(vh); ve_it.is_valid(); ++ve_it)
	{
		each_edge_update_num_count[*ve_it]++;
		Cal_Cost_ForEdge(*ve_it);
	}
}

bool InteractiveViewerWidget::EdgeCollapse(const EdgeCollapseInfo& edge_collapse_info)
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
	else if(mesh.is_collapse_ok(heh_other))
	{
		mesh.set_point(from_p, sp);
		mesh.collapse(heh_other);
		collapse_ok = true;
		v_result = from_p;
	}
	
	if (collapse_ok)
	{
		// 更新 mesh 数据
		// buildIndex();
		// mesh.update_normals();
		Q2v[v_result] = Q_result;
		Update_Local_Variable(v_result);
	}

	return collapse_ok;
}


void InteractiveViewerWidget::NewEdgeCollapse(const EdgeCollapseInfo& edge_collapse_info)
{
	// 最优点
	auto eh = edge_collapse_info.eh;
	auto eV = edge_collapse_info.opt_point;
	std::vector<OpenMesh::HalfedgeHandle> involved_hes;

	auto he = mesh.s_halfedge_handle(eh, 0);
	auto v0 = mesh.from_vertex_handle(he);
	auto v1 = mesh.to_vertex_handle(he);

	// 存储此时要被改变了的边，存储在edge_que内的如果用到这些边了，那就是无效的
	for (auto ve_it = mesh.ve_iter(v0); ve_it.is_valid(); ++ve_it)
	{
		be_changed.insert((*ve_it));
	}
	for (auto ve_it = mesh.ve_iter(v1); ve_it.is_valid(); ++ve_it)
	{
		be_changed.insert((*ve_it));
	}

	for (const auto& vhe : mesh.voh_range(v0)) {
		if (vhe.to() != v1 && vhe.next().to() != v1) {
			involved_hes.push_back(vhe.next());
		}
	}
	for (const auto& vhe : mesh.voh_range(v1)) {
		if (vhe.to() != v0 && vhe.next().to() != v0) {
			involved_hes.push_back(vhe.next());
		}
	}
	mesh.delete_vertex(v0);
	mesh.delete_vertex(v1);
	auto new_v = mesh.add_vertex(eV);
	Q2v[new_v] = Q2v[v0] + Q2v[v1];
	for (const auto& inv : involved_hes) {
		auto f = mesh.add_face(new_v, mesh.from_vertex_handle(inv), mesh.to_vertex_handle(inv));
	}
	
	// buildIndex();
	mesh.update_normals();
	Update_Local_Variable(new_v);
}



// 实现 QEM 简化
void InteractiveViewerWidget::QEMSimplifyMesh()
{
	/*
	SetSimplifyLevel(1);

	
	// 最开始先更新网格的各种 normal，清理中间变量
	mesh.update_normals();
	Q2v.clear();
	swap(empty, edge_que);
	each_edge_update_num_count.clear();

	// 首先计算 mesh 所有点处的 Q 值放入 Q2v 内
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		Cal_Q_ForVertex(*v_it);
	}

	// 计算每条边的 cost
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		each_edge_update_num_count[*e_it] = 0;
		Cal_Cost_ForEdge(*e_it);
	}

	int collapse_num = 0;
	while ((mesh.n_vertices() > simplify_mesh_vertex) && collapse_num < 36)
	{
		// 取出最小的 cost 边
		EdgeCollapseInfo edgeinfo = edge_que.top();
		edge_que.pop();
		if (edgeinfo.cost > cost_min) break;

		auto eh = edgeinfo.eh;
		// 当弹出的元素满足是最新被更新的边才有效
		auto it = std::find(be_changed.begin(), be_changed.end(), eh);
		if (it == be_changed.end())
		{
			NewEdgeCollapse(edgeinfo);
			// std::cout << " "<< eh.idx();
			mesh.update_normals();
			++collapse_num;
			
		}
	}
	Q2v.clear();
	swap(empty, edge_que);
	each_edge_update_num_count.clear();

	mesh.garbage_collection();
	mesh.update_normals();
	*/

	Mesh_Simplification(mesh, 0.3);
	
	mesh_vector.push_back(mesh); mesh_vector_index += 1;
	emit set_edit_undo_enable_viewer_signal(true);
	emit set_edit_redo_enable_viewer_signal(false);
	updateGL();
}

void InteractiveViewerWidget::NewQEMSimplifyMesh()
{
	SetSimplifyLevel(1);

	Q2v.clear();
	swap(empty, edge_que);
	each_edge_update_num_count.clear();

	// 首先计算 mesh 所有点处的 Q 值放入 Q2v 内
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		Cal_Q_ForVertex(*v_it);
	}

	// 计算每条边的 cost
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		each_edge_update_num_count[*e_it] = 0;
		Cal_Cost_ForEdge(*e_it);
	}

	int collapse_num = 0;
	while ((mesh.n_vertices() > simplify_mesh_vertex) && collapse_num < 100)
	{
		// 取出最小的 cost 边
		EdgeCollapseInfo edgeinfo = edge_que.top();
		edge_que.pop();
		auto eh = edgeinfo.eh;
		if (each_edge_update_num_count[eh] == edgeinfo.count)
		{
			EdgeCollapse(edgeinfo);
			++collapse_num;
		}
	}

	Q2v.clear();
	swap(empty, edge_que);
	each_edge_update_num_count.clear();

	mesh.garbage_collection();
	
	mesh_vector.push_back(mesh); mesh_vector_index += 1;
	emit set_edit_undo_enable_viewer_signal(true);
	emit set_edit_redo_enable_viewer_signal(false);
	updateGL();
}


void InteractiveViewerWidget::Mesh_Simplification(Mesh& mesh, float ratio)
{
	assert(ratio >= 0 && ratio <= 1);
	int itr_num = (1.0f - ratio) * mesh.n_vertices();

	// 记录每个点处的 Q 值
	auto Q2v = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, Eigen::Matrix4f>(mesh);
	// 记录每个点从 VertexHandle 到四维向量
	auto v = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, Eigen::Vector4f>(mesh);
	// 配合 vto_flag vfrom_flag 判断点对是否还有效
	auto flag = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, int>(mesh); 
	// p 记录的就是文章中的面的 a,b,c,d
	auto p = OpenMesh::makeTemporaryProperty<OpenMesh::FaceHandle, Eigen::Vector4f>(mesh);

	// 求每个面的表达式
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		Mesh::Normal nor = mesh.normal(*f_it);

		// 代入平面上任意一点，求得 d
		Mesh::Point tp = mesh.point(f_it->halfedge().to());
		double d = -nor | tp;
		p[*f_it][0] = nor[0];
		p[*f_it][1] = nor[1];
		p[*f_it][2] = nor[2];
		p[*f_it][3] = static_cast<float>(d);
	}

	// 根据上面求得的 p 值去求每个点的 Q 值，直接用 p*pT 求和
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		Eigen::Matrix4f Q_tmp;
		Q_tmp.setZero();
		for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(*v_it); vf_it.is_valid(); ++vf_it)
		{
			Q_tmp += p[*vf_it] * (p[*vf_it].transpose());
		}
		Q2v[*v_it] = Q_tmp;
		auto p_v = mesh.point(*v_it);
		v[*v_it][0] = p_v[0];v[*v_it][1] = p_v[1];v[*v_it][2] = p_v[2];
		v[*v_it][3] = 1.0f;
		flag[*v_it] = 0;
	}

	// 使用优先队列每次取出 cost 的最低值
	std::priority_queue <AUX_EdgeCollapse, std::vector<AUX_EdgeCollapse>, 
		std::less<AUX_EdgeCollapse>> edge_queue;
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		// 每条边两侧的 Q 值相加
		Eigen::Matrix4f Q_plus = Q2v[e_it->v0()] + Q2v[e_it->v1()];

		// 求偏导时直接修改最后一行
		Eigen::Matrix4f A = Q_plus;
		Eigen::Vector4f b(0.0f, 0.0f, 0.0f, 1.0f);
		A(3, 0) = 0.0f;
		A(3, 1) = 0.0f;
		A(3, 2) = 0.0f;
		A(3, 3) = 1.0f;

		// 求解出最优位置
		Eigen::Vector4f opt_point = A.colPivHouseholderQr().solve(b);

		AUX_EdgeCollapse edge_collapse;
		edge_collapse.halfeh = e_it->halfedge(0);            // 记录这个边的一条半边
		edge_collapse.cost = opt_point.transpose() * Q_plus * opt_point;  // vT Q v 即为 cost

		Mesh::Point newpt(opt_point[0], opt_point[1], opt_point[2]);
		edge_collapse.newpt = newpt;
		edge_collapse.vto = e_it->halfedge(0).to();
		edge_collapse.vfrom = e_it->halfedge(0).from();
		edge_collapse.Q_new = Q_plus;
		// 按照 cost 进行排序
		edge_queue.push(edge_collapse);
	}


	// 进行迭代
	int count = 0;
	while (count < itr_num)
	{
		// 取出影响最小的边
		AUX_EdgeCollapse min_cost_edge = edge_queue.top();
		edge_queue.pop();

		// 如果这条边的点的句柄已经被删除，跳过
		if (mesh.status(min_cost_edge.vfrom).deleted() || mesh.status(min_cost_edge.vto).deleted())
			continue;
		if (min_cost_edge.vto_flag != flag[min_cost_edge.vto] || 
			min_cost_edge.vfrom_flag != flag[min_cost_edge.vfrom])
			continue;


		Mesh::VertexHandle to_vh;
		if (mesh.is_collapse_ok(min_cost_edge.halfeh))
		{
			// vto 句柄依然存在，向其收缩
			mesh.collapse(min_cost_edge.halfeh);
			to_vh = min_cost_edge.vto;
			// 标记这两个点被合并了一次
			flag[min_cost_edge.vto]++;
			flag[min_cost_edge.vfrom]++;
		}
		// 去收缩对半边
		else if (mesh.is_collapse_ok(mesh.opposite_halfedge_handle(min_cost_edge.halfeh)))
		{
			mesh.collapse(mesh.opposite_halfedge_handle(min_cost_edge.halfeh));
			to_vh = min_cost_edge.vfrom;
			flag[min_cost_edge.vto]++;
			flag[min_cost_edge.vfrom]++;
		}
		else
		{
			continue;
		}

		// 上面并没有收缩到最优点，这里收缩完之后进行一次位置移动至最优点
		mesh.set_point(to_vh, min_cost_edge.newpt);
		Q2v[to_vh] = min_cost_edge.Q_new;
		v[to_vh][0] = min_cost_edge.newpt[0];
		v[to_vh][1] = min_cost_edge.newpt[1];
		v[to_vh][2] = min_cost_edge.newpt[2];
		v[to_vh][3] = 1.0f;

		// 更新这个点附近的所有边的 cost 计算
		for (Mesh::VertexOHalfedgeIter vh_it = mesh.voh_iter(to_vh); vh_it.is_valid(); ++vh_it)
		{
			Mesh::VertexHandle to_v = vh_it->to();
			Eigen::Matrix4f Q_plus = min_cost_edge.Q_new + Q2v[to_v];
			Eigen::Matrix4f A = Q_plus;
			Eigen::Vector4f b(0.0f, 0.0f, 0.0f, 1.0f);
			A(3, 0) = 0.0f;
			A(3, 1) = 0.0f;
			A(3, 2) = 0.0f;
			A(3, 3) = 1.0f;

			// 求解出最优位置
			Eigen::Vector4f opt_point = A.colPivHouseholderQr().solve(b);
			
			AUX_EdgeCollapse edge_collapse;
			edge_collapse.halfeh = *vh_it;
			edge_collapse.cost = opt_point.transpose() * Q_plus * opt_point;

			Mesh::Point newpt(opt_point[0], opt_point[1], opt_point[2]);
			edge_collapse.newpt = newpt;
			edge_collapse.vto = to_v;
			edge_collapse.vto_flag = flag[to_v];
			edge_collapse.vfrom = to_vh;
			edge_collapse.vfrom_flag = flag[to_vh];
			edge_collapse.Q_new = Q_plus;
			edge_queue.push(edge_collapse);
		}
		count++;

	}
	mesh.garbage_collection();
}