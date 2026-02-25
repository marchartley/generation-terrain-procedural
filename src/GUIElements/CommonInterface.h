#ifndef COMMONINTERFACE_H
#define COMMONINTERFACE_H


#include <QWidget>
#include <QMetaObject>
#include <vector>
#include "GUIElements/InterfaceUtils.h"
#include <QPushButton>
#include <QCheckBox>
#include <QRadioButton>
#include <QLineEdit>
#include <QComboBox>
#include <QDial>
#include <QTimer>
#include <optional>
#include <variant>

#include <iostream>

#define DEFINE_SET_ON_FUNCTION(FUNCTION_NAME, WIDGET_TYPE, SIGNAL_NAME) \
template <class Callable, typename... Args> \
    auto FUNCTION_NAME(Callable&& callback, Args&&... args) { \
        addConnection<WIDGET_TYPE>(&WIDGET_TYPE::SIGNAL_NAME, std::forward<Callable>(callback), std::forward<Args>(args)...); \
        return this; \
}

#define DEFINE_SET_ON_SUBWIDGET_FUNCTION(FUNCTION_NAME, WIDGET_TYPE, SUBWIDGET_NAME, SIGNAL_NAME) \
template <class Callable, typename... Args> \
    auto FUNCTION_NAME(Callable&& callback, Args&&... args) { \
        QMetaObject::Connection connection = QObject::connect(SUBWIDGET_NAME, &WIDGET_TYPE::SIGNAL_NAME, \
                           std::forward<Callable>(callback), \
                           std::forward<Args>(args)...); \
        connections.push_back(connection); \
        return this; \
}

class UIElement : public QObject {
    Q_OBJECT
public:
    UIElement(QWidget* widget);
    virtual ~UIElement();

    template <class WidgetType, typename SignalType, typename Functor>
    void addConnection(SignalType signal, Functor&& functor) {
        auto* w = qobject_cast<WidgetType*>(element);
        if (!w) throw std::logic_error("Wrong widget type");
        connections.push_back(QObject::connect(w, signal, this, std::forward<Functor>(functor)));
    }

    UIElement* block() {
        this->blockSignals(true);
        this->element->blockSignals(true);
        return this;
    }
    UIElement* unblock() {
        this->blockSignals(false);
        this->element->blockSignals(false);
        return this;
    }
    /*
    template <class WidgetType, typename SignalType, typename Callable, typename... Args>
    void addConnection(SignalType signal, Callable&& slotFunction, Args&&... args) {
        WidgetType* castedWidget = dynamic_cast<WidgetType*>(element);
        if (castedWidget) {
            QMetaObject::Connection connection = QObject::connect(castedWidget, signal,
                                                                  std::forward<Callable>(slotFunction),
                                                                  std::forward<Args>(args)...);
            connections.push_back(connection);
        } else {
            throw std::out_of_range("There was a problem with a UIElement. It appears that an event has a connection with a wrong QWidget type...");
        }
    }*/

    void setName(const std::string& name);
    const std::string& getName() const;


    virtual QWidget* getWidget() const;
    QWidget* get() const { return getWidget(); }
    void cleanupConnections();

public Q_SLOTS:
    virtual void update();

protected:
    QWidget* element;
    std::vector<QMetaObject::Connection> connections;
    unsigned int ID;
    std::string name;
};


class InterfaceUI : public UIElement {
    Q_OBJECT
public:
    InterfaceUI(QLayout* layout, bool tight = true, std::string title = "");
    ~InterfaceUI();

    QGroupBox* box() const;

    UIElement* add(UIElement* element, std::string name = "");
    std::vector<UIElement*> add(std::vector<UIElement*> elements);
    std::vector<UIElement*> add(std::vector<std::pair<UIElement*, std::string>> elementsAndNames);
    UIElement* add(QLayout* layout, std::string name = "");
    UIElement* find(const std::string& name);
    InterfaceUI *clear();

    std::vector<UIElement*> elements;
    std::vector<std::string> names;
    std::string title;

public Q_SLOTS:
    void update();
};


#include "GUIElements/FancySlider.h"
#include "GUIElements/HierarchicalListWidget.h"
#include "GUIElements/ComboboxElement.h"


class LabelElement : public UIElement {
public:
    LabelElement(const std::string& text);

    QLabel* label();

    LabelElement* setText(const std::string& newText);
    std::string getText();
};

class ButtonElement : public UIElement {
public:
    ButtonElement(const std::string& label);
    ButtonElement(const std::string& label, std::function<void(void)> onClick);

    QPushButton* button();

    ButtonElement* setOnRepeat(std::function<void(void)> onRepeatFunction, int delay_ms = 100);

    DEFINE_SET_ON_FUNCTION(setOnClick, QPushButton, clicked)
    DEFINE_SET_ON_FUNCTION(setOnPressed, QPushButton, pressed)
    DEFINE_SET_ON_FUNCTION(setOnRelease, QPushButton, released)

protected:
    QTimer* pressedTimer = nullptr;
    bool currentRepeatFunctionFinished = true;
};

class SliderElement : public UIElement {
    Q_OBJECT
public:
    SliderElement(const std::string& label, float valMin, float valMax, float multiplier, Qt::Orientation orientation = Qt::Horizontal);
    SliderElement(const std::string& label, float valMin, float valMax, float multiplier, float& binded, Qt::Orientation orientation = Qt::Horizontal);

    SliderElement* setValue(float newValue) { slider()->setfValue(newValue); return this; }
    float value() const { return slider()->getfValue(); }

    FancySlider* slider() const;
    QLabel* label() const;

    DEFINE_SET_ON_SUBWIDGET_FUNCTION(setOnValueChanged, FancySlider, _slider, floatValueChanged)
    /*
    template <class Callable, typename... Args>
    SliderElement* setOnValueChanged(Callable&& callback, Args&&... args) {
        QMetaObject::Connection connection = QObject::connect(_slider, &FancySlider::floatValueChanged,
                                                              std::forward<Callable>(callback),
                                                              std::forward<Args>(args)...);
        connections.push_back(connection);
        return this;
    }
    */

    SliderElement* bindTo(float& value);

public Q_SLOTS:
    void update();

protected:
    QLabel* _label = nullptr;
    FancySlider* _slider = nullptr;
    std::optional<std::reference_wrapper<float>> boundVariable;

    float defaultValue = 0.f;
};


class RangeSliderElement : public UIElement {
    Q_OBJECT
public:
    RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, Qt::Orientation orientation = Qt::Horizontal);
    RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, float& bindedMin, float& bindedMax, Qt::Orientation orientation = Qt::Horizontal);

    RangeSliderElement* setMinValue(float valueMin) { slider()->setMinValue(valueMin); return this; }
    RangeSliderElement* setMaxValue(float valueMax) { slider()->setMaxValue(valueMax); return this; }
    RangeSliderElement* setMinMaxValues(float valueMin, float valueMax) { return this->setMinValue(valueMin)->setMaxValue(valueMax); }

    RangeSlider* slider();
    QLabel* label();

    DEFINE_SET_ON_SUBWIDGET_FUNCTION(setOnValueChanged, RangeSlider, _slider, alt_valueChanged)
    /*
    template <class Callable, typename... Args>
    RangeSliderElement* setOnValueChanged(Callable&& callback, Args&&... args) {
        QMetaObject::Connection connection = QObject::connect(_slider, &RangeSlider::alt_valueChanged,
                                                              std::forward<Callable>(callback),
                                                              std::forward<Args>(args)...);
        connections.push_back(connection);
        return this;
    }
    */

    RangeSliderElement* bindTo(float& valueMin, float& valueMax);

public Q_SLOTS:
    void update();

protected:
    QLabel* _label = nullptr;
    RangeSlider* _slider = nullptr;
    std::optional<std::reference_wrapper<float>> boundVariableMin;
    std::optional<std::reference_wrapper<float>> boundVariableMax;
};

/*
class FloatInputElement : public UIElement {
    Q_OBJECT
public:
    FloatInputElement(const std::string& label);
    FloatInputElement(const std::string& label, float& binded);
    FloatInputElement(const std::string& label, std::function<void(float)> onChange);


};*/

class CheckboxElement : public UIElement {
    Q_OBJECT
public:
    CheckboxElement(const std::string& label);
    CheckboxElement(const std::string& label, bool& binded);
    CheckboxElement(const std::string& label, std::function<void(bool)> onCheck);

    QCheckBox* checkBox();

    CheckboxElement* setChecked(bool checked) { checkBox()->setChecked(checked); return this; }

    DEFINE_SET_ON_FUNCTION(setOnChecked, QCheckBox, toggled)

    CheckboxElement* bindTo(bool& value);

public Q_SLOTS:
    void update();

protected:
    std::optional<std::reference_wrapper<bool>> boundVariable;
};

class RadioButtonElement : public UIElement {
    Q_OBJECT
public:
    RadioButtonElement(const std::string& label);
    RadioButtonElement(const std::string& label, bool& binded);
    RadioButtonElement(const std::string& label, const std::function<void(bool)>& onCheck);

    QRadioButton* radioButton() const;

    RadioButtonElement* setChecked(bool checked) { radioButton()->setChecked(checked); return this; }

    DEFINE_SET_ON_FUNCTION(setOnChecked, QRadioButton, toggled)

    RadioButtonElement* bindTo(bool& value);

    bool checked() const { return radioButton()->isChecked(); }

public Q_SLOTS:
    void update();

protected:
    std::optional<std::reference_wrapper<bool>> boundVariable;
};

class TextEditElement : public UIElement {
    Q_OBJECT
public:
    TextEditElement(const std::string& text, std::string label = "");
    TextEditElement(const std::string& text, std::string label, std::string &binded);

    QLineEdit* lineEdit();
    std::string getText() { return lineEdit()->text().toStdString(); }

    //    DEFINE_SET_ON_FUNCTION(setOnReturnPressed, QLineEdit, returnPressed);

    DEFINE_SET_ON_SUBWIDGET_FUNCTION(setOnReturnPressed, QLineEdit, _lineEdit, returnPressed)
    /*
    template <class Callable, typename... Args>
    TextEditElement* setOnReturnPressed(Callable&& callback, Args&&... args) {
        QMetaObject::Connection connection = QObject::connect(_lineEdit, &QLineEdit::returnPressed,
                                                              std::forward<Callable>(callback),
                                                              std::forward<Args>(args)...);
        connections.push_back(connection);
        return this;
    }
    */
    TextEditElement* setOnTextChange(std::function<void(const std::string&)> func);

    TextEditElement* bindTo(std::string &value);

public Q_SLOTS:
    void update();

protected:
    QLabel* _label = nullptr;
    QLineEdit* _lineEdit = nullptr;
    std::optional<std::reference_wrapper<std::string>> boundVariable;
};


class AngleElement : public UIElement {
    Q_OBJECT
public:
    AngleElement(const std::string& label = "");
    AngleElement(const std::string& label, float &binded);

    QDial* dial() const;
    float getAngle() const { return this->value(); }
    float value() const { return (float) dial()->value(); }

    DEFINE_SET_ON_SUBWIDGET_FUNCTION(setOnValueChanged, QDial, _dial, valueChanged)

    AngleElement* bindTo(float& value);

public Q_SLOTS:
    void update();

protected:
    QLabel* _label = nullptr;
    QDial* _dial = nullptr;
    std::optional<std::reference_wrapper<float>> boundVariable;
};


#include "GUIElements/qtcolorpicker.h"
#include "DataStructure/Vector3.h"
class ColorPickerElement : public UIElement {
    Q_OBJECT
public:
    ColorPickerElement(const std::string& label);
    ColorPickerElement(const std::string& label, Vector3& currentSelection);

    QtColorPicker* colorPicker() const;

    ColorPickerElement* setOnSelectionChanged(std::function<void(const Vector3&)> func);

    ColorPickerElement* bindTo(Vector3& indexSelected);

    Vector3 getSelection() const;

public Q_SLOTS:
    void update();

    //protected:
public:
    std::optional<std::reference_wrapper<Vector3>> boundColor;

    QLabel* _label;
    QtColorPicker* _colorPicker;

protected:
    static Vector3 qColorToVec3(const QColor& col);
    static QColor vec3ToQColor(const Vector3& col);
};

class HierarchicalListUI : public UIElement {
public:
    HierarchicalListUI();

    HierarchicalListWidget* hierarchicalList();

    HierarchicalListUI* setSelectionMode(QAbstractItemView::SelectionMode mode);

    DEFINE_SET_ON_FUNCTION(setOnItemSelectionChanged, HierarchicalListWidget, itemSelectionChanged);

    HierarchicalListUI* clear();
    HierarchicalListUI* addItem(HierarchicalListWidgetItemBase* item);
    HierarchicalListUI* setCurrentItems(std::vector<int> selectedIDs);
    std::vector<HierarchicalListWidgetItemBase*> selectedItems();

    template <class T>
    HierarchicalListUI* removeItem(const T& itemToRemove)
    {
        this->hierarchicalList()->removeItem(itemToRemove);
        return this;
    }
};


InterfaceUI* createHorizontalGroupUI(std::vector<UIElement*> widgets);
InterfaceUI* createVerticalGroupUI(std::vector<UIElement*> widgets);
InterfaceUI* createMultiColumnGroupUI(std::vector<UIElement*> widgets, int nbColumns = 2);

#undef DEFINE_SET_ON_FUNCTION


#endif // COMMONINTERFACE_H
