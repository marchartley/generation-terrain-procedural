#include "InterfaceUtils.h"

#include <QLabel>
#include <QBoxLayout>
#include <QCheckBox>
#include <QApplication>

#include <unordered_map>

/*
QGroupBox* createSliderGroup(const std::string& label, QSlider* slider, bool makeItSmall)
{
    QLabel* lab = new QLabel(QString::fromStdString(label));
    QBoxLayout* layout;
    if (makeItSmall)
        layout = new QHBoxLayout;
    else
        layout = new QVBoxLayout;
    layout->setMargin(0);
    QGroupBox* group = new QGroupBox;
    layout->addWidget(lab);
    layout->addWidget(slider);
    group->setLayout(layout);

    return group;
}
QGroupBox* createMultipleSliderGroup(std::vector<std::pair<std::string, QSlider*>> labelsAndSliders)
{
    QGridLayout* layout = new QGridLayout();
    QGroupBox* group = new QGroupBox;
    int row = 0;
    for (const auto& labAndSlid : labelsAndSliders) {
        QLabel* lab = new QLabel(QString::fromStdString(std::get<0>(labAndSlid)));
        layout->addWidget(lab, row, 0);
        layout->addWidget(std::get<1>(labAndSlid), row, 1);
        row++;
    }
    group->setLayout(layout);

    return group;
}
QGroupBox* createVerticalGroup(std::vector<QWidget*> widgets)
{
    QVBoxLayout* layout = new QVBoxLayout;
    QGroupBox* group = new QGroupBox;
    for (QWidget*& w : widgets)
        layout->addWidget(w);
    group->setLayout(layout);
    return group;
}

QGroupBox* createOptionalSlider(RangeSlider *slider, std::string checkboxLabel, bool activatedByDefault, std::function<void(bool, RangeSlider*)> onToggleCallback)
{
}

QGroupBox *createMultipleSliderGroupWithCheckbox(std::vector<std::tuple<std::string, QSlider*, QCheckBox*>> labelsAndSlidersAndActivables)
{
    QGridLayout* layout = new QGridLayout();
    QGroupBox* group = new QGroupBox;
    int row = 0;
    for (const auto& labAndSlidAndAct : labelsAndSlidersAndActivables) {
        QLabel* lab = new QLabel(QString::fromStdString(std::get<0>(labAndSlidAndAct)));
        layout->addWidget(lab, row, 0);
        layout->addWidget(std::get<1>(labAndSlidAndAct), row, 1);
        layout->addWidget(std::get<2>(labAndSlidAndAct), row, 2);
//        layout->addWidget(std::get<0>(std::get<1>(labAndSlidAndAct)), row, 1);
//        layout->addWidget(std::get<1>(std::get<1>(labAndSlidAndAct)), row, 2);
        row++;
    }
    group->setLayout(layout);
    group->resize(group->size());
    return group;
}

QGroupBox *createHorizontalGroup(std::vector<QWidget *> widgets)
{
    QHBoxLayout* layout = new QHBoxLayout;
    QGroupBox* group = new QGroupBox;
    for (QWidget*& w : widgets)
        layout->addWidget(w);
    group->setLayout(layout);
    return group;
}

QGroupBox *createMultiColumnGroup(std::vector<QWidget *> widgets, int nbColumns)
{
    QGridLayout* layout = new QGridLayout;
    QGroupBox* group = new QGroupBox;
//    for (QWidget*& w : widgets)
    for (size_t i = 0; i < widgets.size(); i++) {
        auto& w = widgets[i];
        int row = int (i / nbColumns);
        int col = i % nbColumns;
        layout->addWidget(w, row, col);
    }
    group->setLayout(layout);
    return group;
}
*/
