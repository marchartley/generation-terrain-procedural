#ifndef COMMONINTERFACE_H
#define COMMONINTERFACE_H


#include <QWidget>
#include <QLayout>
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
#include <QDoubleSpinBox>
#include <optional>
#include <variant>

#include <iostream>

#define DEFINE_SET_ON_FUNCTION(FUNCTION_NAME, WIDGET_TYPE, SIGNAL_NAME) \
template <class Callable, typename... Args> \
    auto& FUNCTION_NAME(Callable&& callback, Args&&... args) { \
        addConnection<WIDGET_TYPE>(&WIDGET_TYPE::SIGNAL_NAME, std::forward<Callable>(callback), std::forward<Args>(args)...); \
        return *this; \
}

#define DEFINE_SET_ON_SUBWIDGET_FUNCTION(FUNCTION_NAME, WIDGET_TYPE, SUBWIDGET_NAME, SIGNAL_NAME) \
template <class Callable, typename... Args> \
    auto& FUNCTION_NAME(Callable&& callback, Args&&... args) { \
        QMetaObject::Connection connection = QObject::connect(SUBWIDGET_NAME, &WIDGET_TYPE::SIGNAL_NAME, \
                           std::forward<Callable>(callback), \
                           std::forward<Args>(args)...); \
        connections.push_back(connection); \
        return *this; \
}

class UIElement : public QWidget {
    Q_OBJECT
public:
    UIElement(QWidget* widget = nullptr);
    virtual ~UIElement();

    template <class WidgetType, typename SignalType, typename Functor>
    void addConnection(SignalType signal, Functor&& functor) {
        auto* w = qobject_cast<WidgetType*>(element);
        if (!w) throw std::logic_error("Wrong widget type");
        connections.push_back(QObject::connect(w, signal, this, std::forward<Functor>(functor)));
    }

    UIElement& block() {
        this->blockSignals(true);
        this->element->blockSignals(true);
        return *this;
    }
    UIElement& unblock() {
        this->blockSignals(false);
        this->element->blockSignals(false);
        return *this;
    }

    UIElement& setName(const std::string& name);
    const std::string& getName() const;


    virtual QWidget* getWidget() const;
    QWidget* get() const { return getWidget(); }
    void cleanupConnections();


    enum LAYOUT { HORIZONTAL, VERTICAL, GRID };

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
    InterfaceUI(LAYOUT layout = VERTICAL, const std::string& title = "");
    InterfaceUI(LAYOUT layout, bool tight, const std::string& title = "");
    ~InterfaceUI();

    QGroupBox* box() const;

    template<class T, class... Args>
    T& add(Args&&... args)
    {
        static_assert(std::is_base_of_v<UIElement, T>,
                      "T must derive from UIElement");

        auto element = new T(std::forward<Args>(args)...);

        T* raw = element;
        this->add(raw);

        return *raw;
    }
    UIElement* add(UIElement* element, const std::string& name = "");
    void add(const std::vector<UIElement*>& elements);
    InterfaceUI& clear();

    InterfaceUI* setTight(bool tight);

    std::vector<UIElement*> elements;
    std::vector<std::string> names;
    std::string title;

public Q_SLOTS:
    void update();
};

class LabelElement : public UIElement {
public:
    LabelElement(const std::string& text);

    QLabel* label();

    LabelElement& setText(const std::string& newText);
    std::string getText();
};

#include "GUIElements/FancySlider.h"
#include "GUIElements/HierarchicalListWidget.h"
#include "GUIElements/ComboboxElement.h"



class ButtonElement : public UIElement {
public:
    ButtonElement(const std::string& label);
    ButtonElement(const std::string& label, std::function<void(void)> onClick);

    QPushButton* button();

    ButtonElement& setOnRepeat(std::function<void(void)> onRepeatFunction, int delay_ms = 100, int startup_time = 1000);

    DEFINE_SET_ON_FUNCTION(setOnClick, QPushButton, clicked)
    DEFINE_SET_ON_FUNCTION(setOnPressed, QPushButton, pressed)
    DEFINE_SET_ON_FUNCTION(setOnRelease, QPushButton, released)

protected:
    QTimer* pressedTimer = nullptr;
    bool currentRepeatFunctionFinished = true;
};


class SliderElement : public InterfaceUI {
    Q_OBJECT
public:
    SliderElement(const std::string& label, float valMin, float valMax, float multiplier, UIElement::LAYOUT orientation = HORIZONTAL);
    SliderElement(const std::string& label, float valMin, float valMax, float multiplier, float& binded, UIElement::LAYOUT orientation = HORIZONTAL);
    SliderElement(const std::string& label, float valMin, float valMax, float multiplier, std::function<void(float)> onValueChanged, UIElement::LAYOUT orientation = HORIZONTAL);
    SliderElement(const std::string& label, float valMin, float valMax, float multiplier, float& binded, std::function<void(float)> onValueChanged, UIElement::LAYOUT orientation = HORIZONTAL);

    SliderElement& setValue(float newValue) { slider()->setfValue(newValue); return *this; }
    float value() const { return slider()->getfValue(); }
    float getValue() const { return value(); }

    FancySlider* slider() const;
    LabelElement* label() const;

    DEFINE_SET_ON_SUBWIDGET_FUNCTION(setOnValueChanged, FancySlider, _slider, floatValueChanged)

    SliderElement& bindTo(float& value);

public Q_SLOTS:
    void update();

protected:
    LabelElement* _label = nullptr;
    FancySlider* _slider = nullptr;
    std::optional<std::reference_wrapper<float>> boundVariable;

    float defaultValue = 0.f;
};


class RangeSliderElement : public InterfaceUI {
    Q_OBJECT
public:
    RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, UIElement::LAYOUT orientation = HORIZONTAL);
    RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, float& bindedMin, float& bindedMax, UIElement::LAYOUT orientation = HORIZONTAL);
    RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, std::function<void(float, float)> onChange, UIElement::LAYOUT orientation = HORIZONTAL);
    RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, float& bindedMin, float& bindedMax, std::function<void(float, float)> onChange, UIElement::LAYOUT orientation = HORIZONTAL);

    RangeSliderElement& setMinValue(float valueMin) { slider()->setMinValue(valueMin); return *this; }
    RangeSliderElement& setMaxValue(float valueMax) { slider()->setMaxValue(valueMax); return *this; }
    RangeSliderElement& setMinMaxValues(float valueMin, float valueMax) { return this->setMinValue(valueMin).setMaxValue(valueMax); }

    RangeSlider* slider();
    LabelElement* label();

    DEFINE_SET_ON_SUBWIDGET_FUNCTION(setOnValueChanged, RangeSlider, _slider, alt_valueChanged)

    RangeSliderElement& bindTo(float& valueMin, float& valueMax);

public Q_SLOTS:
    void update();

protected:
    LabelElement* _label = nullptr;
    RangeSlider* _slider = nullptr;
    std::optional<std::reference_wrapper<float>> boundVariableMin;
    std::optional<std::reference_wrapper<float>> boundVariableMax;
};


class CheckboxElement : public UIElement {
    Q_OBJECT
public:
    CheckboxElement(const std::string& label);
    CheckboxElement(const std::string& label, bool& binded);
    CheckboxElement(const std::string& label, std::function<void(bool)> onCheck);
    CheckboxElement(const std::string& label, bool& binded, std::function<void(bool)> onCheck);

    QCheckBox* checkBox();

    CheckboxElement& setChecked(bool checked) { checkBox()->setChecked(checked); return *this; }

    DEFINE_SET_ON_FUNCTION(setOnChecked, QCheckBox, toggled)

    CheckboxElement& bindTo(bool& value);

public Q_SLOTS:
    void update();

protected:
    std::optional<std::reference_wrapper<bool>> boundVariable;
};

class TextEditElement : public InterfaceUI {
    Q_OBJECT
public:
    TextEditElement(const std::string& text, std::string label = "");
    TextEditElement(const std::string& text, std::string label, std::string &binded);

    QLineEdit* lineEdit();
    std::string getText() { return lineEdit()->text().toStdString(); }
    TextEditElement& setText(const std::string& newText);

    DEFINE_SET_ON_SUBWIDGET_FUNCTION(setOnReturnPressed, QLineEdit, _lineEdit, returnPressed)
    TextEditElement& setOnTextChange(std::function<void(const std::string&)> func);

    TextEditElement& bindTo(std::string &value);

public Q_SLOTS:
    void update();

protected:
    LabelElement* _label = nullptr;
    QLineEdit* _lineEdit = nullptr;
    std::optional<std::reference_wrapper<std::string>> boundVariable;
};


class FloatInputElement : public InterfaceUI {
    Q_OBJECT
public:
    FloatInputElement(const std::string& label = "");
    FloatInputElement(const std::string& label, float& binded);
    FloatInputElement(const std::string& label, std::function<void(float)> onChange);

    QDoubleSpinBox* spinbox() const { return this->_spinbox; }
    float getValue() const { return spinbox()->value(); }
    FloatInputElement& setValue(float newValue);

    // DEFINE_SET_ON_SUBWIDGET_FUNCTION(setOnValueChanged, QDoubleSpinBox, _spinbox, valueChanged)
    auto& setOnValueChanged(std::function<void(const float&)> func); // QOverload<double>::of(&QDoubleSpinBox::valueChanged)
    template <class Callable, typename... Args> \
    FloatInputElement& setOnValueChanged(Callable&& callback, Args&&... args) {
        QMetaObject::Connection connection = QObject::connect(_spinbox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                           std::forward<Callable>(callback),
                           std::forward<Args>(args)...);
        connections.push_back(connection);
        return *this;
    }

    FloatInputElement& bindTo(float& value);

public Q_SLOTS:
    void update();

protected:
    LabelElement* _label = nullptr;
    QDoubleSpinBox* _spinbox = nullptr;
    std::optional<std::reference_wrapper<float>> boundVariable;


};


class AngleElement : public InterfaceUI {
    Q_OBJECT
public:
    AngleElement(const std::string& label = "");
    AngleElement(const std::string& label, float &binded);

    QDial* dial() const;
    float getAngle() const { return this->value(); }
    float value() const { return (float) dial()->value(); }
    AngleElement& setAngle(float newAngle);

    DEFINE_SET_ON_SUBWIDGET_FUNCTION(setOnValueChanged, QDial, _dial, valueChanged)

    AngleElement& bindTo(float& value);

public Q_SLOTS:
    void update();

protected:
    LabelElement* _label = nullptr;
    QDial* _dial = nullptr;
    std::optional<std::reference_wrapper<float>> boundVariable;
};


#include "GUIElements/qtcolorpicker.h"
#include "DataStructure/Vector3.h"
class ColorPickerElement : public InterfaceUI {
    Q_OBJECT
public:
    ColorPickerElement(const std::string& label);
    ColorPickerElement(const std::string& label, Vector3& currentSelection);

    QtColorPicker* colorPicker() const;

    ColorPickerElement& setOnSelectionChanged(std::function<void(const Vector3&)> func);

    ColorPickerElement& bindTo(Vector3& indexSelected);

    Vector3 getSelection() const;

public Q_SLOTS:
    void update();

    //protected:
public:
    std::optional<std::reference_wrapper<Vector3>> boundColor;

    LabelElement* _label;
    QtColorPicker* _colorPicker;

protected:
    static Vector3 qColorToVec3(const QColor& col);
    static QColor vec3ToQColor(const Vector3& col);
};

class HierarchicalListUI : public UIElement {
public:
    HierarchicalListUI();

    HierarchicalListWidget* hierarchicalList();

    HierarchicalListUI& setSelectionMode(QAbstractItemView::SelectionMode mode);

    DEFINE_SET_ON_FUNCTION(setOnItemSelectionChanged, HierarchicalListWidget, itemSelectionChanged);

    HierarchicalListUI& clear();
    HierarchicalListUI& addItem(HierarchicalListWidgetItemBase* item);
    HierarchicalListUI& setCurrentItems(std::vector<int> selectedIDs);
    std::vector<HierarchicalListWidgetItemBase*> selectedItems();

    template <class T>
    HierarchicalListUI& removeItem(const T& itemToRemove)
    {
        this->hierarchicalList()->removeItem(itemToRemove);
        return *this;
    }
};


InterfaceUI* createHorizontalGroup(std::vector<UIElement*> widgets, bool tight = true);
InterfaceUI* createVerticalGroup(std::vector<UIElement*> widgets, bool tight = true);
InterfaceUI* createMultiColumnGroup(std::vector<UIElement*> widgets, int nbColumns = 2, bool tight = true);








template <class T>
class RadioButtonElement : public UIElement {
    // Q_OBJECT
public:
    // RadioButtonElement(const std::string& label);
    RadioButtonElement(const std::string& label, T associatedValue, T& binded);
    RadioButtonElement(const std::string& label, const std::function<void(bool)>& onCheck);

    QRadioButton* radioButton() const;

    RadioButtonElement& setChecked(bool checked) { radioButton()->setChecked(checked); return *this; }

    DEFINE_SET_ON_FUNCTION(setOnChecked, QRadioButton, toggled)

    RadioButtonElement& bindTo(T& value);

    bool checked() const { return radioButton()->isChecked(); }

// public Q_SLOTS:
    void update();

protected:
    std::optional<std::reference_wrapper<T>> boundVariable;
    std::optional<T> associatedValue;
};

template <class T>
RadioButtonElement<T>::RadioButtonElement(const std::string& label, T associatedValue, T& binded)
    : UIElement(new QRadioButton(QString::fromStdString(label)))
{
    this->associatedValue = associatedValue;
    bindTo(binded);
}
template <class T>
RadioButtonElement<T>::RadioButtonElement(const std::string& label, const std::function<void(bool)>& onCheckedFunction)
    : UIElement(new QRadioButton(QString::fromStdString(label)))
{
    this->setOnChecked(onCheckedFunction);
}

template <class T>
QRadioButton* RadioButtonElement<T>::radioButton() const
{
    return qobject_cast<QRadioButton*>(getWidget());
}

template <class T>
RadioButtonElement<T>& RadioButtonElement<T>::bindTo(T &value)
{
    radioButton()->setChecked(value == associatedValue.value());
    boundVariable = value;
    setOnChecked([this](bool checked) {
        if (checked && boundVariable) {
            boundVariable->get() = this->associatedValue.value();
        }
    });
    return *this;
}

template <class T>
void RadioButtonElement<T>::update()
{
    if (this->boundVariable.has_value())
        radioButton()->setChecked(this->boundVariable->get() == associatedValue.value());
    return UIElement::update();
}

#undef DEFINE_SET_ON_FUNCTION


#endif // COMMONINTERFACE_H
