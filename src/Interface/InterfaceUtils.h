#ifndef INTERFACEUTILS_H
#define INTERFACEUTILS_H

#include <string>
#include <QSlider>
#include <QGroupBox>

QGroupBox* createSliderGroup(std::string label, QSlider* slider, bool makeItSmall = false);
QGroupBox* createMultipleSliderGroup(std::map<std::string, QSlider*> labelsAndSliders);
QGroupBox* createVerticalGroup(std::vector<QWidget*> widgets);

#endif // INTERFACEUTILS_H
