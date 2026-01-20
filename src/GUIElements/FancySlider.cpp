#include "FancySlider.h"

#include <cmath>
#include <QStyleOption>
#include <QToolTip>
#include <QPainter>
#include <QDebug>

#include "Utils/Utils.h"

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
    // Give space for ticks labels
    setMinimumHeight(QSlider::sizeHint().height() * 2);
    setMinimumWidth(QSlider::sizeHint().width() * 2);
}
FancySlider::~FancySlider()
{

}

float FancySlider::getfValue()
{
    return this->value() * this->multiplier;
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

void FancySlider::addTicks(const std::vector<std::pair<float, std::string> > &posAndLabels)
{
    this->ticksPositionsAndLabels = posAndLabels;
}

void FancySlider::paintEvent(QPaintEvent *ev)
{
    QSlider::paintEvent(ev);
    if (!this->ticksPositionsAndLabels.empty()) {

        QPainter painter = QPainter(this);
        painter.setPen(QPen(Qt::black));

        auto font_metrics = QFontMetrics(this->font());


        auto rect = this->geometry();
        if (this->orientation() == Qt::Horizontal) {
            for (const std::pair<float, std::string>& tick : this->ticksPositionsAndLabels) {
                float pos = tick.first;
                QString lab = QString::fromStdString(tick.second);

                auto font_width = font_metrics.boundingRect(lab).width();
                auto font_height = font_metrics.boundingRect(lab).height();

                auto x_pos = rect.width() * interpolation::linear(pos / multiplier, this->minimum(), this->maximum());
                auto y_pos = rect.height() * 1.f;

                auto label_x_pos = std::max(0.f, std::min(x_pos - (font_width * .5f), float(rect.width() - font_width)));
                painter.drawLine(QPoint(x_pos, rect.height() * .25f), QPoint(x_pos, rect.height() * .5f));
                painter.drawText(QPoint(label_x_pos, y_pos), lab);
            }

        } else if (this->orientation() == Qt::Vertical) {
            /*auto vertical_x_pos = rect.width() - font_width - 5;
            auto vertical_y_pos = rect.height() * 0.75;
            painter.drawText(QPoint(rect.width() / 2.0 - font_width / 2.0, rect.height() - 5),
                              lab);*/
        } else {

        }

        painter.drawRect(rect);
    }
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
