#include "RangeSlider.h"

#include <QMouseEvent>
#include "qmimedata.h"
#include "qdrag.h"
#include "qwidgetaction.h"
#include "qapplication.h"
#include "qpixmap.h"
#include "qcursor.h"
#include "qguiapplication.h"
#include "qdir.h"
#include <QProxyStyle>

#include <iostream>

class SliderProxy : public QProxyStyle
{
public:
  int pixelMetric ( PixelMetric metric, const QStyleOption * option = 0, const QWidget * widget = 0 ) const
  {
    switch(metric) {
    case PM_SliderThickness  : return 25;
    case PM_SliderLength     : return 25;
    default                  : return (QProxyStyle::pixelMetric(metric,option,widget));
    }
  }
};

RangeSlider::RangeSlider(Qt::Orientation orientation, float min, float max, float multiplier, QWidget *parent)
  : FancySlider(orientation, min, max, multiplier, parent)
{
  //styling
  setAcceptDrops(true);
  SliderProxy *aSliderProxy = new SliderProxy();

  //hard coded path to image :/ sorry
//  QString path = QDir::fromNativeSeparators(this->pathToHandleImage);
  setStyleSheet("QSlider::handle { image: none; }"); //image: url(" + path + "); }");
  setStyle(aSliderProxy);

  //setting up the alternate handle
  min_handle = new RangeSliderHandle(0.f, this);
  addAction(new QWidgetAction(min_handle));
//  min_handle->move(this->pos().x() + this->width()- min_handle->width(), this->pos().y() );

  //setting up the alternate handle
  max_handle = new RangeSliderHandle(1.f, this);
  addAction(new QWidgetAction(max_handle));
//  max_handle->move(this->pos().x() + max_handle->width(), this->pos().y() );

}

RangeSliderHandle::RangeSliderHandle(float val, RangeSlider *_parent)
  : QLabel(_parent), val(val)
{
  parent = _parent;
  handleActivated = false;
  filter = new SuperEventFilter(parent);

  //styling
  setAcceptDrops(true);
  //hard coded path to image :/ sorry
  QPixmap pix = QPixmap::fromImage(QImage(parent->pathToHandleImage));
  pix =  pix.scaled(25, 25, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
  setPixmap(pix);
//  this->updatePos();
}

void RangeSlider::sliderChange(SliderChange change)
{

}

void RangeSlider::resizeEvent(QResizeEvent *e)
{
    FancySlider::resizeEvent(e);
    min_handle->updatePos();
    max_handle->updatePos();
}

float RangeSlider::min_value()
{
  return min_handle->value();
}
float RangeSlider::max_value()
{
  return max_handle->value();
}

void RangeSlider::setMinValue(float value)
{
  min_handle->setValue(value);
}
void RangeSlider::setMaxValue(float value)
{
  max_handle->setValue(value);
}

void RangeSlider::mouseReleaseEvent(QMouseEvent *mouseEvent)
{
  if (mouseEvent->button() == Qt::LeftButton)
  {
      min_handle->show();
      min_handle->handleActivated = false;

      max_handle->show();
      max_handle->handleActivated = false;
  }
//  mouseEvent->accept();
}

void RangeSlider::mousePressEvent(QMouseEvent *mouseEvent)
{
  if (mouseEvent->button() == Qt::LeftButton)
  {
      QPoint posCursor = this->mapFromGlobal(QCursor::pos());
      QPoint point_clicked(posCursor.x(), 0);
     if (std::abs(min_handle->pos().x() - point_clicked.x()) < std::abs(max_handle->pos().x() - point_clicked.x()))
     {
         min_handle->handleActivated = true;
     } else {
         max_handle->handleActivated = true;
     }
  }
  //  mouseEvent->accept();
}

void RangeSlider::mouseMoveEvent(QMouseEvent *event)
{
    this->alt_update();
}

float RangeSlider::getValueUnderCursor(QPoint cur)
{
    double width = this->width(), position = cur.x();
    double value = position/width;
    return value;
//    double range = this->maximum() - this->minimum();
//    std::cout << value << std::endl;
//    return this->minimum() + (value * range);
}

void RangeSlider::alt_update()
{
   QPoint posCursor = this->mapFromGlobal(QCursor::pos());
   QPoint point_clicked(posCursor.x(), 0);
  if (min_handle->handleActivated)
  {
      min_handle->setValue(this->getValueUnderCursor(posCursor));
  } else if (max_handle->handleActivated) {
      max_handle->setValue(this->getValueUnderCursor(posCursor));
  }
    this->constraintMinMaxValues();
  Q_EMIT alt_valueChanged(min_value(), max_value());
  Q_EMIT minValueChanged(min_value());
  Q_EMIT maxValueChanged(max_value());
}

void RangeSlider::constraintMinMaxValues()
{
    if (min_handle->val > max_handle->val)
    {
        float tmp_val = min_handle->val;
        bool tmp_hold = min_handle->handleActivated;
        min_handle->val = max_handle->val;
        min_handle->handleActivated = max_handle->handleActivated;
        max_handle->val = tmp_val;
        max_handle->handleActivated = tmp_hold;
    }
}
/*
void RangeSliderHandle::mousePressEvent(QMouseEvent *mouseEvent)
{
    parent->mousePressEvent(mouseEvent);
//    std::cout << "Pressed child" << std::endl;
  qGuiApp->installEventFilter(filter);
  parent->clearFocus();
}*/

bool SuperEventFilter::eventFilter(QObject* obj, QEvent* event)
{
   /* if (obj != grandParent)
        return false;*/
    std::cout << event->type() << std::endl;
  switch(event->type())
  {
  case QEvent::MouseButtonRelease:
      grandParent->mouseReleaseEvent((QMouseEvent*)event);
    qGuiApp->removeEventFilter(this);
//    return true;
    break;
  case QEvent::MouseButtonPress:
      grandParent->mousePressEvent((QMouseEvent*)event);
//      break;
  case QEvent::MouseMove:
    grandParent->alt_update();
//    return true;
    break;
  default:
    return QObject::eventFilter(obj, event);
  }
  return false;
}

void RangeSliderHandle::setValue(float value)
{
    this->val = std::clamp(value, 0.f, 1.f);
    this->updatePos();
    /**/
}

float RangeSliderHandle::value()
{
    return ((this->val - this->parent->minimum()) / (this->parent->maximum() - this->parent->minimum()) + parent->minimum()) / this->parent->multiplier;
    /*
  double width = parent->width(), position = pos().x();
  double value = position/width;
  double range = parent->maximum() - parent->minimum();
  return parent->minimum() + (value * range);*/
}

void RangeSliderHandle::updatePos()
{
    this->move(this->val * parent->width() - this->width() / 2, 0);/*
    double width = parent->width(); //, position = pos().x();
    double range = parent->maximum() - parent->minimum();
    float location = (this->val - parent->minimum())/range;
    location = location * width;
    move(y(),location);*/
}
void RangeSlider::Reset()
{
  int horBuffer = (min_handle->width());
  QPoint myPos = mapToGlobal(pos());
  QPoint point(myPos.x() + width() - horBuffer, myPos.y()- min_handle->height());
  point = min_handle->mapFromParent(point);

  min_handle->move(point);
  min_handle->show();
  min_handle->raise();

}
