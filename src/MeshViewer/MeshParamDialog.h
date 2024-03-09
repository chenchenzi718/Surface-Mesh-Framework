#ifndef MESHPROCESSING_MESHPARAMDIALOG_H
#define MESHPROCESSING_MESHPARAMDIALOG_H

#include <QDialog>
#include <QtGui>
#include <QtWidgets>


/*
*  创建了整个框架的左侧 tab QP 以及按钮 Print Mesh Information
*/

class MeshParamDialog : public QDialog
{
	Q_OBJECT
public:
	MeshParamDialog(QWidget* parent=0);
	~MeshParamDialog();

	QSize sizeHint()
	{
		QRect rect = QApplication::desktop()->screenGeometry();
		return QSize( int( rect.width()*0.15), rect.height() );
	}

private:
	QTabWidget* tabWidget;

signals:
	void print_info_signal();

private:
	QWidget* Basic_Operation_And_Information;
	QScrollArea *view_BOI;

	QLabel* leftLabel_BOI;
	QPushButton* print_info;

private:
	void create_Basic_Operation_Information_Widget();

private:
	void initDialog();
	void createWidget();
	void createLayout();

};

#endif