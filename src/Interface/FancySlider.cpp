#include "FancySlider.h"

#include <cmath>
#include <QStyleOption>
#include <QToolTip>

FancySlider::FancySlider(float multiplier, QWidget *parent)
        : QSlider(parent), multiplier(multiplier)
{
    connect(this, SIGNAL(valueChanged(int)),
                this, SLOT(notifyValueChanged(int)));
}
FancySlider::FancySlider(Qt::Orientation orientation, float min, float max, float multiplier, QWidget *parent)
        : QSlider(orientation, parent), multiplier(multiplier)
{
    connect(this, SIGNAL(valueChanged(int)),
                this, SLOT(notifyValueChanged(int)));
    this->setRange(std::round(min / multiplier), std::round(max / multiplier));
}
FancySlider::~FancySlider()
{

}
void FancySlider::setfValue(float val)
{
    int a = int(std::round(val / multiplier));
    this->setValue(a);
}
void FancySlider::setfRange(float min, float max)
{
    this->setRange(int(std::round(min / multiplier)), int(std::round(max / multiplier)));
}
void FancySlider::notifyValueChanged(int value) {
    double doubleValue = value * multiplier;
    Q_EMIT floatValueChanged(doubleValue);
}
void FancySlider::sliderChange(SliderChange change)
{
    QSlider::sliderChange(change);

    if (change == QAbstractSlider::SliderValueChange )
    {
        QStyleOptionSlider opt;
        initStyleOption(&opt);

        QRect sr = style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderHandle, this);
        QPoint bottomRightCorner = sr.bottomLeft();

        QToolTip::showText(mapToGlobal( QPoint( bottomRightCorner.x(), bottomRightCorner.y() ) ), QString::number((double)value() * multiplier), this);
    }
}
