#ifndef INTERFACEUTILS_H
#define INTERFACEUTILS_H

#include <string>
#include <QSlider>
#include <QGroupBox>
#include <QCheckBox>
#include <functional>
#include "GUIElements/RangeSlider.h"

QGroupBox* createSliderGroup(const std::string& label, QSlider* slider, bool makeItSmall = false);
QGroupBox* createMultipleSliderGroup(std::vector<std::pair<std::string, QSlider *> > labelsAndSliders);
QGroupBox* createMultipleSliderGroupWithCheckbox(std::vector<std::tuple<std::string, QSlider*, QCheckBox*>> labelsAndSlidersAndActivables);
QGroupBox* createVerticalGroup(std::vector<QWidget*> widgets);
QGroupBox* createHorizontalGroup(std::vector<QWidget*> widgets);
QGroupBox* createMultiColumnGroup(std::vector<QWidget*> widgets, int nbColumns = 2);
QGroupBox* createOptionalSlider(RangeSlider* slider, std::string checkboxLabel = "Activer", bool activatedByDefault = true,
                                std::function<void(bool, RangeSlider*)> onToggleCallback = []([[maybe_unused]] bool active, [[maybe_unused]] RangeSlider* slider) -> void {});

//QWidget* stickTo(QWidget* widget, QWidget* container, float x, float y, float w, float h, bool useAbsolutePosition = true);
//void clearLayout(QLayout* layout, bool deleteWidgets = true);

#endif // INTERFACEUTILS_H
