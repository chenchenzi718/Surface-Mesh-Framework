#include"Aux_QEMWidget.h"


AUX_QEMWidget::AUX_QEMWidget(QWidget* parent)
    :QWidget(parent)
{
    // �����������
    sliderTitle = new QLabel("ratio of simplification", this);

    // ��������
    slider = new QSlider(Qt::Horizontal);
    slider->setRange(0, 100); // ����Ϊ0-100�ķ�Χ
    slider->setValue(50);

    // ������ǩ������ʾ����ֵ
    sliderValueLabel = new QLabel("0.5", this); // ��ʼֵ�뻬��һ��

    // ���ӻ����valueChanged�źŵ��ۺ����Ը��±�ǩ
    connect(slider, &QSlider::valueChanged, this, &AUX_QEMWidget::label_value_change);
    connect(slider, &QSlider::valueChanged, this, &AUX_QEMWidget::qem_slider_value_changed);

    // ����һ�������������ֵ��ǩ�Ĳ���
    QHBoxLayout* sliderLayout = new QHBoxLayout();
    sliderLayout->addWidget(slider);
    sliderLayout->addWidget(sliderValueLabel);

    // ����һ����ֱ���֣�������ӱ��⣬Ȼ����ӻ��鲼��
    QVBoxLayout* titledSliderLayout = new QVBoxLayout();
    titledSliderLayout->addWidget(sliderTitle);
    titledSliderLayout->addLayout(sliderLayout);

    // ������
    this->setLayout(titledSliderLayout);
}


void AUX_QEMWidget::label_value_change(int value)
{
    double mappedValue = value / 100.0; // ����Χ��0-100ӳ�䵽0-1
    sliderValueLabel->setText(QString::number(mappedValue, 'f', 2));
}


void AUX_QEMWidget::qem_slider_value_changed(int value)
{
    emit qem_slider_value_changed_signal(value);
}