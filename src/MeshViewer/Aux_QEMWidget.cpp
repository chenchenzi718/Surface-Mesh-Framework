#include"Aux_QEMWidget.h"


AUX_QEMWidget::AUX_QEMWidget(QWidget* parent)
    :QWidget(parent)
{
    // 创建滑块标题
    sliderTitle = new QLabel("ratio of simplification", this);

    // 创建滑块
    slider = new QSlider(Qt::Horizontal);
    slider->setRange(0, 100); // 设置为0-100的范围
    slider->setValue(50);

    // 创建标签用于显示滑块值
    sliderValueLabel = new QLabel("0.5", this); // 初始值与滑块一致

    // 连接滑块的valueChanged信号到槽函数以更新标签
    connect(slider, &QSlider::valueChanged, this, &AUX_QEMWidget::label_value_change);
    connect(slider, &QSlider::valueChanged, this, &AUX_QEMWidget::qem_slider_value_changed);

    // 创建一个包含滑块和其值标签的布局
    QHBoxLayout* sliderLayout = new QHBoxLayout();
    sliderLayout->addWidget(slider);
    sliderLayout->addWidget(sliderValueLabel);

    // 创建一个垂直布局，首先添加标题，然后添加滑块布局
    QVBoxLayout* titledSliderLayout = new QVBoxLayout();
    titledSliderLayout->addWidget(sliderTitle);
    titledSliderLayout->addLayout(sliderLayout);

    // 主布局
    this->setLayout(titledSliderLayout);
}


void AUX_QEMWidget::label_value_change(int value)
{
    double mappedValue = value / 100.0; // 将范围从0-100映射到0-1
    sliderValueLabel->setText(QString::number(mappedValue, 'f', 2));
}


void AUX_QEMWidget::qem_slider_value_changed(int value)
{
    emit qem_slider_value_changed_signal(value);
}