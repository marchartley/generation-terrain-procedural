#include "InterfaceUtils.h"

#include <QLabel>
#include <QBoxLayout>
#include <QCheckBox>
#include <QApplication>

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
/*
QWidget *stickTo(QWidget* widget, QWidget *container, float x, float y, float w, float h, bool useAbsolutePosition)
{
    widget->setParent(container);
    widget->installEventFilter(container);
//    widget->eventFilter()
}
void clearLayout(QLayout* layout, bool deleteWidgets)
{
    while (QLayoutItem* item = layout->takeAt(0))
    {
        if (deleteWidgets)
        {
            if (QWidget* widget = item->widget())
                widget->deleteLater();
        }
        if (QLayout* childLayout = item->layout())
            clearLayout(childLayout, deleteWidgets);
        delete item;
    }
}
*/

QGroupBox* createOptionalSlider(RangeSlider *slider, std::string checkboxLabel, bool activatedByDefault, std::function<void(bool, RangeSlider*)> onToggleCallback)
{
    /*
    QHBoxLayout* layout = new QHBoxLayout();
    QGroupBox* group = new QGroupBox;

    QCheckBox* checkbox = new QCheckBox(QString::fromStdString(checkboxLabel));
    checkbox->setChecked(activatedByDefault);
    QObject::connect(checkbox, &QCheckBox::toggled,
                     slider, [&](bool activation) { return onToggleCallback(activation, slider); });
//    QObject::connect(checkbox, &QCheckBox::toggled, slider, [&](bool active) {
//        if (active) {
//            WidgetActivationEvent event;
//            QApplication::postEvent(slider, &event);
//        } else {
//            WidgetDesactivationEvent event;
//            QApplication::postEvent(slider, &event);
//        }
//    });

    layout->addWidget(slider);
    layout->addWidget(checkbox);

    group->setLayout(layout);

    return group;
    */
}

QGroupBox *createMultipleSliderGroupWithCheckbox(std::map<std::string, std::pair<QSlider *, QCheckBox *>> labelsAndSlidersAndActivables)
{
    QGridLayout* layout = new QGridLayout();
    QGroupBox* group = new QGroupBox;
    int row = 0;
    for (const auto& labAndSlidAndAct : labelsAndSlidersAndActivables) {
        QLabel* lab = new QLabel(QString::fromStdString(std::get<0>(labAndSlidAndAct)));
        layout->addWidget(lab, row, 0);
        layout->addWidget(std::get<0>(std::get<1>(labAndSlidAndAct)), row, 1);
        layout->addWidget(std::get<1>(std::get<1>(labAndSlidAndAct)), row, 2);
        row++;
    }
    group->setLayout(layout);
    group->resize(group->size());
    return group;
}
