#ifndef INTERFACEUTILS_H
#define INTERFACEUTILS_H

#include <string>
#include <QSlider>
#include <QGroupBox>
#include <QCheckBox>
#include <functional>
#include "Interface/RangeSlider.h"

QGroupBox* createSliderGroup(std::string label, QSlider* slider, bool makeItSmall = false);
QGroupBox* createMultipleSliderGroup(std::map<std::string, QSlider*> labelsAndSliders);
QGroupBox* createMultipleSliderGroupWithCheckbox(std::map<std::string, std::pair<QSlider*, QCheckBox*>> labelsAndSlidersAndActivables);
QGroupBox* createVerticalGroup(std::vector<QWidget*> widgets);
QGroupBox* createOptionalSlider(RangeSlider* slider, std::string checkboxLabel = "Activer", bool activatedByDefault = true,
                                std::function<void(bool, RangeSlider*)> onToggleCallback = [](bool active, RangeSlider* slider) -> void {});

//QWidget* stickTo(QWidget* widget, QWidget* container, float x, float y, float w, float h, bool useAbsolutePosition = true);
//void clearLayout(QLayout* layout, bool deleteWidgets = true);

#endif // INTERFACEUTILS_H
