#ifndef COMMONINTERFACE_H
#define COMMONINTERFACE_H

#include <QWidget>
#include <QMetaObject>
#include <vector>
#include "Interface/InterfaceUtils.h"
#include "Interface/FancySlider.h"
#include "Interface/HierarchicalListWidget.h"
#include <QPushButton>
#include <QCheckBox>
#include <QRadioButton>
#include <QLineEdit>
#include <QComboBox>
#include <optional>
#include <variant>

class UIElement;



#define DEFINE_SET_ON_FUNCTION(FUNCTION_NAME, WIDGET_TYPE, SIGNAL_NAME) \
    template <typename Callable, typename... Args> \
    void FUNCTION_NAME(Callable&& callback, Args&&... args) { \
        addConnection<WIDGET_TYPE>(&WIDGET_TYPE::SIGNAL_NAME, std::forward<Callable>(callback), std::forward<Args>(args)...); \
    }

class UIElement : public QObject {
    Q_OBJECT
public:
    UIElement(QWidget* widget);
    ~UIElement();

    template <typename WidgetType, typename SignalType, typename Callable, typename... Args>
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
    }

    void setName(std::string name);
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

class LabelElement : public UIElement {
public:
    LabelElement(std::string text);

    QLabel* label();

    LabelElement* setText(std::string newText);
    std::string getText();
};

class ButtonElement : public UIElement {
public:
    ButtonElement(std::string label);
    ButtonElement(std::string label, std::function<void(void)> onClick);

    QPushButton* button();

    DEFINE_SET_ON_FUNCTION(setOnClick, QPushButton, clicked)
    DEFINE_SET_ON_FUNCTION(setOnPressed, QPushButton, pressed)
    DEFINE_SET_ON_FUNCTION(setOnRelease, QPushButton, released)


};

class SliderElement : public UIElement {
    Q_OBJECT
public:
    SliderElement(std::string label, float valMin, float valMax, float multiplier, Qt::Orientation orientation = Qt::Horizontal);
    SliderElement(std::string label, float valMin, float valMax, float multiplier, float& binded, Qt::Orientation orientation = Qt::Horizontal);

    FancySlider* slider();
    QLabel* label();

    template <typename Callable, typename... Args>
    void setOnValueChanged(Callable&& callback, Args&&... args) {
        QMetaObject::Connection connection = QObject::connect(_slider, &FancySlider::floatValueChanged,
                                                              std::forward<Callable>(callback),
                                                              std::forward<Args>(args)...);
        connections.push_back(connection);
    }

    void bindTo(float& value);

public Q_SLOTS:
    void update();

protected:
    QLabel* _label = nullptr;
    FancySlider* _slider = nullptr;
    std::optional<std::reference_wrapper<float>> boundVariable;
};


class RangeSliderElement : public UIElement {
    Q_OBJECT
public:
    RangeSliderElement(std::string label, float valMin, float valMax, float multiplier, Qt::Orientation orientation = Qt::Horizontal);
    RangeSliderElement(std::string label, float valMin, float valMax, float multiplier, float& bindedMin, float& bindedMax, Qt::Orientation orientation = Qt::Horizontal);

    RangeSlider* slider();
    QLabel* label();

    template <typename Callable, typename... Args>
    void setOnValueChanged(Callable&& callback, Args&&... args) {
        QMetaObject::Connection connection = QObject::connect(_slider, &RangeSlider::alt_valueChanged,
                                                              std::forward<Callable>(callback),
                                                              std::forward<Args>(args)...);
        connections.push_back(connection);
    }

    void bindTo(float& valueMin, float& valueMax);

public Q_SLOTS:
    void update();

protected:
    QLabel* _label = nullptr;
    RangeSlider* _slider = nullptr;
    std::optional<std::reference_wrapper<float>> boundVariableMin;
    std::optional<std::reference_wrapper<float>> boundVariableMax;
};

class CheckboxElement : public UIElement {
    Q_OBJECT
public:
    CheckboxElement(std::string label);
    CheckboxElement(std::string label, bool& binded);
    CheckboxElement(std::string label, std::function<void(bool)> onCheck);

    QCheckBox* checkBox();

    void setChecked(bool checked) { checkBox()->setChecked(checked); }

    DEFINE_SET_ON_FUNCTION(setOnChecked, QCheckBox, toggled)

    void bindTo(bool& value);

public Q_SLOTS:
    void update();

protected:
    std::optional<std::reference_wrapper<bool>> boundVariable;
};

class RadioButtonElement : public UIElement {
    Q_OBJECT
public:
    RadioButtonElement(std::string label);
    RadioButtonElement(std::string label, bool& binded);

    QRadioButton* radioButton();

    DEFINE_SET_ON_FUNCTION(setOnChecked, QRadioButton, toggled)

    void bindTo(bool& value);

public Q_SLOTS:
    void update();

protected:
    std::optional<std::reference_wrapper<bool>> boundVariable;
};

class TextEditElement : public UIElement {
    Q_OBJECT
public:
    TextEditElement(std::string text, std::string label = "");
    TextEditElement(std::string text, std::string label, std::string &binded);

    QLineEdit* lineEdit();
    std::string getText() { return lineEdit()->text().toStdString(); }

//    DEFINE_SET_ON_FUNCTION(setOnReturnPressed, QLineEdit, returnPressed);

    template <typename Callable, typename... Args>
    void setOnReturnPressed(Callable&& callback, Args&&... args) {
        QMetaObject::Connection connection = QObject::connect(_lineEdit, &QLineEdit::returnPressed,
                                                              std::forward<Callable>(callback),
                                                              std::forward<Args>(args)...);
        connections.push_back(connection);
    }
    void setOnTextChange(std::function<void(std::string)> func);

    void bindTo(std::string& value);

public Q_SLOTS:
    void update();

protected:
    QLabel* _label = nullptr;
    QLineEdit* _lineEdit = nullptr;
    std::optional<std::reference_wrapper<std::string>> boundVariable;
};


struct ComboboxLineElement {
//    ComboboxLineElement(std::string label) : label(label) {}
//    ComboboxLineElement(std::string label, int value) : label(label), value(value) {}
//    ComboboxLineElement(std::string label) : label(label) {}
    std::string label;
    int value;
    std::string iconPath;
    std::string otherParameters;
};

class ComboboxElement : public UIElement {
    Q_OBJECT
public:
    ComboboxElement(std::string label);
    ComboboxElement(std::string label, std::vector<ComboboxLineElement> choices);
    ComboboxElement(std::string label, std::vector<ComboboxLineElement> choices, int& currentSelection);

    QComboBox* combobox() const;

    void setOnSelectionChanged(std::function<void(int)> func);

    void bindTo(int& indexSelected);

    ComboboxLineElement getSelection() const;

public Q_SLOTS:
    void update();

//protected:
public:
    std::vector<ComboboxLineElement> choices;
    std::optional<std::reference_wrapper<int>> boundIndex;
//    std::optional<std::reference_wrapper<std::vector<std::string>>> boundValues;

    QLabel* _label;
    QComboBox* _combobox;

//    bool itemsAreImages = false;
};

class HierarchicalListUI : public UIElement {
public:
    HierarchicalListUI();

    HierarchicalListWidget* hierarchicalList();

    HierarchicalListUI* setSelectionMode(QAbstractItemView::SelectionMode mode);

    DEFINE_SET_ON_FUNCTION(setOnItemSelectionChanged, HierarchicalListWidget, itemSelectionChanged);

    HierarchicalListUI* clear();
    HierarchicalListUI* addItem(HierarchicalListWidgetItem* item);
    HierarchicalListUI* setCurrentItems(std::vector<int> selectedIDs);
    QList<QListWidgetItem *> selectedItems();
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
    UIElement* find(std::string name);

    std::vector<UIElement*> elements;
    std::vector<std::string> names;
    std::string title;

public Q_SLOTS:
    void update();
};

InterfaceUI* createHorizontalGroupUI(std::vector<UIElement*> widgets);
InterfaceUI* createVerticalGroupUI(std::vector<UIElement*> widgets);
InterfaceUI* createMultiColumnGroupUI(std::vector<UIElement*> widgets, int nbColumns = 2);

#undef DEFINE_SET_ON_FUNCTION

#endif // COMMONINTERFACE_H
