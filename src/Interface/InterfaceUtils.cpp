#include "InterfaceUtils.h"

#include <QLabel>
#include <QBoxLayout>

QGroupBox* createSliderGroup(std::string label, QSlider* slider, bool makeItSmall)
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
QGroupBox* createMultipleSliderGroup(std::map<std::string, QSlider*> labelsAndSliders)
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
