#ifndef RANGESLIDER_H
#define RANGESLIDER_H
#pragma once

#include "qslider.h"
#include "qlabel.h"

#include "Interface/FancySlider.h"


/*
*  Super sick nasty awesome double handled slider!
*   https://stackoverflow.com/a/26025796/9863298
*   @author Steve
*/
class RangeSliderHandle;

class RangeSlider: public FancySlider
{
  Q_OBJECT
public:
  /** Constructor */
  RangeSlider(Qt::Orientation orientation, float min = 0.f, float max = 1.f, float multiplier = 0.1f, QWidget *parent = 0);

  /** Store the alternate handle for this slider*/
  RangeSliderHandle *min_handle;
  RangeSliderHandle *max_handle;

  /** Overridden mouse release event */
  void mouseReleaseEvent(QMouseEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent* event);

  float getValueUnderCursor(QPoint cur);

  /** Returns the slider value for the alternate handle */
  float min_value();
  float max_value();

  /** Convenience function for setting the value of the alternate handle */
  void setMinValue(float value);
  void setMaxValue(float value);

  void sliderChange(SliderChange change);
  void resizeEvent(QResizeEvent* e);

  /** Resets the alternate handle to the right side of the slider */
  void Reset();

  /** Used to update the position of the alternate handle through the use of an event filter */
  void alt_update();

  void constraintMinMaxValues();

  QString pathToHandleImage = ":/src/assets/handle.png";

Q_SIGNALS:
  /** Constructor */
  void alt_valueChanged(float, float);
  void minValueChanged(float);
  void maxValueChanged(float);
};

class SuperEventFilter : public QObject
{
public:
  /** Constructor */
  SuperEventFilter(RangeSlider *_grandParent)
  {grandParent = _grandParent;};

protected:
  /*
  * overloaded functions for object that inherit from this class
  */
  bool eventFilter(QObject* obj, QEvent* event);

private:
  /** Store the SuperSlider that this event filter is used within. */
  RangeSlider *grandParent;
};

class RangeSliderHandle: public QLabel
{
  Q_OBJECT
public:
  /** Constructor */
  RangeSliderHandle(float val = 0.f, RangeSlider *parent = 0);

  /** An overloaded mousePressevent so that we can start grabbing the cursor and using it's position for the value */
//  void mousePressEvent(QMouseEvent *event);
//  void mouseReleaseEvent(QMouseEvent *event);
//  void mouseMoveEvent(QMouseEvent* event);

  /** Returns the value of this handle with respect to the slider */
  float value();

  /** Maps mouse coordinates to slider values */
  float mapValue();

  void updatePos();

  /** Store the parent as a slider so that you don't have to keep casting it  */
  RangeSlider *parent;

  /** Store a bool to determine if the alternate handle has been activated  */
  bool handleActivated;
  float val;

private:
  /** Store the filter for installation on the qguiapp */
  SuperEventFilter *filter;

public Q_SLOTS:
    /** Sets the value of the handle with respect to the slider */
    void setValue(float value);
};

#endif // RANGESLIDER_H
