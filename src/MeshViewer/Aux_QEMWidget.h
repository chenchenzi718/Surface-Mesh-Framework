#ifndef AUX_QEMWIDGET_D
#define AUX_QEMWIDGET_D

#include <QWidget>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QSlider>
#include <QLabel>

#include"QEMSimplification.h"


// Îª QEM Ìí¼Ó»¬¿é
class AUX_QEMWidget :public QWidget
{
    Q_OBJECT

public:
    AUX_QEMWidget(QWidget* parent = nullptr);

signals:
    void qem_slider_value_changed_signal(int value);

private slots:
    void label_value_change(int value);
    void qem_slider_value_changed(int value);

private:
    QSlider* slider;

    QLabel* sliderTitle, *sliderValueLabel;

};
#endif // !AUX_QEMWIDGET