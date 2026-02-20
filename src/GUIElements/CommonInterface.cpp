#include "CommonInterface.h"

#include <iostream>
#include <QBoxLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>

UIElement::UIElement(QWidget* widget)
    : element(widget)
{

}

QWidget* UIElement::getWidget() const {
    return element;
}

void UIElement::cleanupConnections() {
    for (const auto& connection : connections) {
        QObject::disconnect(connection);
    }
    connections.clear();
}

void UIElement::update()
{
    return;
}


UIElement::~UIElement() {
    cleanupConnections();
    element->deleteLater();
}

void UIElement::setName(std::string name)
{
    this->name = name;
}
const std::string &UIElement::getName() const
{
    return name;
}


LabelElement::LabelElement(std::string text)
    : UIElement(new QLabel(QString::fromStdString(text)))
{}

QLabel* LabelElement::label()
{
    return qobject_cast<QLabel*>(getWidget());
}

LabelElement* LabelElement::setText(std::string newText)
{
    this->label()->setText(QString::fromStdString(newText));
    return this;
}

std::string LabelElement::getText()
{
    return this->label()->text().toStdString();
}



ButtonElement::ButtonElement(const std::string& label)
    : UIElement(new QPushButton(QString::fromStdString(label))) {

    QObject::connect(this->button(), &QObject::destroyed, this, [this](){
        if (pressedTimer) pressedTimer->stop();
    });
}

ButtonElement::ButtonElement(const std::string& label, std::function<void ()> onClick)
    : UIElement(new QPushButton(QString::fromStdString(label)))
{
    this->setOnClick(onClick);
    QObject::connect(this->button(), &QObject::destroyed, this, [this](){
        if (pressedTimer) pressedTimer->stop();
    });
}

QPushButton* ButtonElement::button() {
    return qobject_cast<QPushButton*>(getWidget());
}

ButtonElement* ButtonElement::setOnRepeat(std::function<void ()> onRepeatFunction, int delay_ms)
{
    if (this->pressedTimer != nullptr) delete this->pressedTimer;
    this->pressedTimer = new QTimer(this);
    this->pressedTimer->setInterval(delay_ms);
    QObject::connect(this->pressedTimer, &QTimer::timeout, this, [=](){
        if (this->currentRepeatFunctionFinished) {
            this->currentRepeatFunctionFinished = false;
            onRepeatFunction();
            this->currentRepeatFunctionFinished = true;
        }
        if(this->button() && this->button()->isDown())
            this->pressedTimer->setInterval(delay_ms);
    });
    this->setOnPressed([=]() {
        this->pressedTimer->start(1000.f);
    });
    this->setOnRelease([=]() {
        this->pressedTimer->stop();
    });
    return this;
}




SliderElement::SliderElement(const std::string& label, float valMin, float valMax, float multiplier, Qt::Orientation orientation)
    : UIElement(new QGroupBox)
{
    this->_slider = new FancySlider(orientation, valMin, valMax, multiplier);
    this->_label = new QLabel(QString::fromStdString(label));

    QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty())
        layout->addWidget(_label);
    layout->addWidget(_slider);
    getWidget()->setLayout(layout);

    defaultValue = valMin;

    QObject::connect(this->slider(), &FancySlider::doubleClicked, this, [=]() { this->setValue(defaultValue); });
}

SliderElement::SliderElement(const std::string& label, float valMin, float valMax, float multiplier, float &binded, Qt::Orientation orientation)
    : SliderElement(label, valMin, valMax, multiplier, orientation)
{
    bindTo(binded);
    this->defaultValue = binded;
}

FancySlider* SliderElement::slider() const {
    return this->_slider;
}

QLabel* SliderElement::label() const
{
    return this->_label;
}

SliderElement* SliderElement::bindTo(float &value) {
    _slider->setfValue(value);
    boundVariable = value;
    setOnValueChanged([this](float newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
    return this;
}

void SliderElement::update()
{
    if (this->boundVariable.has_value())
        slider()->setfValue(*this->boundVariable);
}



CheckboxElement::CheckboxElement(const std::string& label)
    : UIElement(new QCheckBox(QString::fromStdString(label)))
{

}

CheckboxElement::CheckboxElement(const std::string& label, bool &binded)
    : CheckboxElement(label)
{
    bindTo(binded);
}

CheckboxElement::CheckboxElement(const std::string& label, std::function<void (bool)> onCheck)
    : CheckboxElement(label)
{
    this->setOnChecked(onCheck);
}

QCheckBox* CheckboxElement::checkBox() {
    return qobject_cast<QCheckBox*>(getWidget());
}

CheckboxElement* CheckboxElement::bindTo(bool &value)
{
    checkBox()->setChecked(value);
    boundVariable = value;
    setOnChecked([this](bool newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
    return this;
}

void CheckboxElement::update()
{
    if (this->boundVariable.has_value())
        checkBox()->setChecked(*this->boundVariable);
}

InterfaceUI::InterfaceUI(QLayout* layout, bool tight, std::string title)
    : UIElement(new QGroupBox), title(title)
{
    if (tight) {
        layout->setSpacing(0);
        layout->setContentsMargins(0, 0, 0, 0);
        get()->setProperty("class", "tight");
    }
    box()->setStyleSheet("QGroupBox.tight{ border: none; }");
    getWidget()->setLayout(layout);
}

InterfaceUI::~InterfaceUI()
{
    /*for (auto& child : elements) {
        child->deleteLater();
    }*/
    elements.clear();
}

QGroupBox* InterfaceUI::box() const
{
    return qobject_cast<QGroupBox*>(getWidget());
}

UIElement* InterfaceUI::add(UIElement* element, std::string name)
{
    box()->layout()->addWidget(element->get());
    element->setName(name);
    this->elements.push_back(element);
    this->names.push_back(name);
    return element;
}

std::vector<UIElement*> InterfaceUI::add(std::vector<UIElement* > elements)
{
    for (auto& element : elements)
        this->add(element);
    return elements;
}

std::vector<UIElement*> InterfaceUI::add(std::vector<std::pair<UIElement* , std::string> > elementsAndNames)
{
    std::vector<UIElement*> justElements(elementsAndNames.size());
    for (size_t i = 0; i < elementsAndNames.size(); i++) {
        auto& element = elementsAndNames[i].first;
        auto& name = elementsAndNames[i].second;
        this->add(element, name);
        justElements[i] = element;
    }
    return justElements;
}

UIElement* InterfaceUI::add(QLayout* layout, std::string name)
{
    InterfaceUI* interface = new InterfaceUI(layout, true, name);
    return this->add(interface, name);
}

UIElement* InterfaceUI::find(std::string name)
{
    for (size_t i = 0; i < names.size(); i++)
        if (names[i] == name)
            return elements[i];
    return nullptr;
}

InterfaceUI* InterfaceUI::clear()
{
    for (auto& child : this->elements) {// box()->children()) {
        child->deleteLater();
    }
    this->elements.clear();
    this->names.clear();
    return this;
}

void InterfaceUI::update()
{
    for (auto& child : this->elements) {// box()->children()) {
        if (auto asUIElement = dynamic_cast<UIElement*>(child)) {
            asUIElement->update();
        }
    }
}

RadioButtonElement::RadioButtonElement(const std::string& label)
    : UIElement(new QRadioButton(QString::fromStdString(label)))
{

}

RadioButtonElement::RadioButtonElement(const std::string& label, bool &binded)
    : RadioButtonElement(label)
{
    bindTo(binded);
}

RadioButtonElement::RadioButtonElement(const std::string& label, const std::function<void (bool)> &onCheck)
    : RadioButtonElement(label)
{
    this->setOnChecked(onCheck);
}


QRadioButton* RadioButtonElement::radioButton() const
{
    return qobject_cast<QRadioButton*>(getWidget());
}

RadioButtonElement* RadioButtonElement::bindTo(bool &value)
{
    radioButton()->setChecked(value);
    boundVariable = value;
    setOnChecked([this](bool newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
    return this;
}

void RadioButtonElement::update()
{
    if (this->boundVariable.has_value())
        radioButton()->setChecked(*this->boundVariable);
}

InterfaceUI* createHorizontalGroupUI(std::vector<UIElement*> widgets)
{
    InterfaceUI* interface = new InterfaceUI(new QHBoxLayout);
    for (UIElement*& w : widgets) {
        interface->add(w);
    }
    return interface;
}

InterfaceUI* createVerticalGroupUI(std::vector<UIElement*> widgets)
{
    InterfaceUI* interface = new InterfaceUI(new QVBoxLayout);
    for (UIElement*& w : widgets) {
        interface->add(w);
    }
    return interface;
}

InterfaceUI* createMultiColumnGroupUI(std::vector<UIElement*> widgets, int nbColumns)
{
    QGridLayout* layout = new QGridLayout;
    // QGroupBox* group = new QGroupBox;
    //    for (QWidget*& w : widgets)
    for (size_t i = 0; i < widgets.size(); i++) {
        auto& w = widgets[i];
        int row = int (i / nbColumns);
        int col = i % nbColumns;

        layout->addWidget(w->get(), row, col);
        // layout->addWidget(w->get(), row, col);
    }
    // group->setLayout(layout);
    return new InterfaceUI(layout); // group;
}

TextEditElement::TextEditElement(std::string text, std::string label)
    : UIElement(new QGroupBox)
{
    this->_lineEdit = new QLineEdit(QString::fromStdString(text));
    this->_label = new QLabel(QString::fromStdString(label));

    QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty())
        layout->addWidget(_label);
    layout->addWidget(_lineEdit);
    getWidget()->setLayout(layout);
}

TextEditElement::TextEditElement(std::string text, std::string label, std::string &binded)
    : TextEditElement(text, label)
{
    this->bindTo(binded);
}

QLineEdit* TextEditElement::lineEdit()
{
    return dynamic_cast<QLineEdit*>(_lineEdit);
}

TextEditElement* TextEditElement::setOnTextChange(std::function<void (std::string)> func)
{
//    this->addConnection()
    QObject::connect(_lineEdit, &QLineEdit::textChanged, this, [=](QString newText){ // /!\ Capture function by value
        func(newText.toStdString());
    });
    return this;
}

TextEditElement* TextEditElement::bindTo(std::string &value)
{
    lineEdit()->setText(QString::fromStdString(value));
    boundVariable = value;
    setOnTextChange([this](std::string newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
    return this;
}

void TextEditElement::update()
{
    lineEdit()->setText(QString::fromStdString(*this->boundVariable));
}


AngleElement::AngleElement(const std::string& label)
    : UIElement(new QGroupBox)
{
    this->_dial = new QDial();
    this->_dial->setWrapping(true); this->_dial->setMinimum(0); this->_dial->setMaximum(360);
    this->_label = new QLabel(QString::fromStdString(label));

    QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty())
        layout->addWidget(_label);
    layout->addWidget(_dial);
    getWidget()->setLayout(layout);
}

AngleElement::AngleElement(const std::string& label, float &binded)
    : AngleElement(label)
{
    this->bindTo(binded);
}

QDial *AngleElement::dial() const
{
    return this->_dial;
}

AngleElement *AngleElement::bindTo(float &value)
{
    dial()->setValue((int)value);
    boundVariable = value;
    setOnValueChanged([this](float newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
    return this;
}

void AngleElement::update()
{
    this->dial()->setValue(*boundVariable);
}


RangeSliderElement::RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, Qt::Orientation orientation)
    : UIElement(new QGroupBox)
{
    this->_slider = new RangeSlider(orientation, valMin, valMax, multiplier);
    this->_label = new QLabel(QString::fromStdString(label));

    QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty()) {
        layout->addWidget(_label);
    }
    layout->addWidget(_slider);
    getWidget()->setLayout(layout);
}

RangeSliderElement::RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, float &bindedMin, float &bindedMax, Qt::Orientation orientation)
    : RangeSliderElement(label, valMin, valMax, multiplier, orientation)
{
    this->bindTo(bindedMin, bindedMax);
}

RangeSlider* RangeSliderElement::slider()
{
    return _slider;
}

QLabel* RangeSliderElement::label()
{
    return _label;
}

RangeSliderElement* RangeSliderElement::bindTo(float &valueMin, float &valueMax)
{
    slider()->setMinValue(valueMin);
    slider()->setMaxValue(valueMax);
    boundVariableMin = valueMin;
    boundVariableMax = valueMax;
    setOnValueChanged([this](float newMin, float newMax) {
        if (boundVariableMin) {
            boundVariableMin->get() = newMin;
        }
        if (boundVariableMax) {
            boundVariableMax->get() = newMax;
        }
    });
    return this;
}

void RangeSliderElement::update()
{
    slider()->setMinValue(*this->boundVariableMin);
    slider()->setMaxValue(*this->boundVariableMax);
}



ColorPickerElement::ColorPickerElement(const std::string& label)
    : UIElement(new QGroupBox)
{
    this->_colorPicker = new QtColorPicker(this->get());
    this->_label = new QLabel(QString::fromStdString(label));

    QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty())
        layout->addWidget(_label);
    layout->addWidget(_colorPicker);
    getWidget()->setLayout(layout);
}

ColorPickerElement::ColorPickerElement(const std::string& label, Vector3 &currentSelection)
    : ColorPickerElement(label)
{
    this->bindTo(currentSelection);
}

QtColorPicker* ColorPickerElement::colorPicker() const
{
    return _colorPicker;
}

ColorPickerElement* ColorPickerElement::setOnSelectionChanged(std::function<void (const Vector3&)> func)
{
    QObject::connect(_colorPicker, &QtColorPicker::colorChanged, this, [=](QColor color) {
        func(ColorPickerElement::qColorToVec3(color));
    });
    return this;
}

ColorPickerElement* ColorPickerElement::bindTo(Vector3& colorSelected)
{
    boundColor = colorSelected;
    colorPicker()->setCurrentColor(ColorPickerElement::vec3ToQColor(colorSelected));
    this->setOnSelectionChanged([&](Vector3 color) {
        this->boundColor->get() = color;
    });
    return this;
}

Vector3 ColorPickerElement::getSelection() const
{
    return ColorPickerElement::qColorToVec3(this->colorPicker()->currentColor());
}

void ColorPickerElement::update()
{
    colorPicker()->setCurrentColor(ColorPickerElement::vec3ToQColor(*this->boundColor));
}

Vector3 ColorPickerElement::qColorToVec3(const QColor &col)
{
    return Vector3(col.redF(), col.greenF(), col.blueF());
}

QColor ColorPickerElement::vec3ToQColor(const Vector3 &col)
{
    return QColor(col.r(), col.g(), col.b());
}


HierarchicalListUI::HierarchicalListUI()
    : UIElement(new HierarchicalListWidget())
{}

HierarchicalListWidget *HierarchicalListUI::hierarchicalList()
{
    return qobject_cast<HierarchicalListWidget*>(this->getWidget());
}

HierarchicalListUI* HierarchicalListUI::setSelectionMode(QAbstractItemView::SelectionMode mode)
{
    this->hierarchicalList()->setSelectionMode(mode);
    return this;
}


HierarchicalListUI* HierarchicalListUI::clear()
{
    this->hierarchicalList()->clear();
    return this;
}

HierarchicalListUI* HierarchicalListUI::addItem(HierarchicalListWidgetItemBase* item)
{
    this->hierarchicalList()->addItem(item);
    return this;
}

HierarchicalListUI* HierarchicalListUI::setCurrentItems(std::vector<int> selectedIDs)
{
    this->hierarchicalList()->setCurrentItems(selectedIDs);
    return this;
}

std::vector<HierarchicalListWidgetItemBase *> HierarchicalListUI::selectedItems()
{
    auto qItems = this->hierarchicalList()->selectedItems();
    std::vector<HierarchicalListWidgetItemBase*> items;

    for (int i = 0; i < qItems.size(); i++) {
        items.push_back(dynamic_cast<HierarchicalListWidgetItemBase*>(qItems[i]));
    }
    return items;
}
