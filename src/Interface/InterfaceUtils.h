#ifndef INTERFACEUTILS_H
#define INTERFACEUTILS_H

#include <string>
#include <QSlider>
#include <QGroupBox>

QGroupBox* createSliderGroup(std::string label, QSlider* slider, bool makeItSmall = false);
QGroupBox* createMultipleSliderGroup(std::map<std::string, QSlider*> labelsAndSliders);
QGroupBox* createVerticalGroup(std::vector<QWidget*> widgets);

//QWidget* stickTo(QWidget* widget, QWidget* container, float x, float y, float w, float h, bool useAbsolutePosition = true);
//void clearLayout(QLayout* layout, bool deleteWidgets = true);

#endif // INTERFACEUTILS_H
