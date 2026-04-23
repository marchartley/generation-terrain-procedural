#include "CommonInterface.h"

#include <iostream>
#include <QBoxLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>

UIElement::UIElement(QWidget* widget)
    : element(widget)
{
    // element->setParent(this);
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
    element->update();
    return this->get()->update();
}



UIElement::~UIElement() {
    cleanupConnections();
    element->deleteLater();
}

UIElement& UIElement::setName(const std::string& name)
{
    this->name = name;
    return *this;
}
const std::string& UIElement::getName() const
{
    return name;
}


LabelElement::LabelElement(const std::string& text)
    : UIElement(new QLabel(QString::fromStdString(text)))
{}

QLabel* LabelElement::label()
{
    return qobject_cast<QLabel*>(getWidget());
}

LabelElement& LabelElement::setText(const std::string& newText)
{
    this->label()->setText(QString::fromStdString(newText));
    return *this;
}

std::string LabelElement::getText()
{
    return this->label()->text().toStdString();
}



ButtonElement::ButtonElement(const std::string& label)
    : UIElement(new QPushButton(QString::fromStdString(label)))
{
    this->pressedTimer = new QTimer(this);
    QObject::connect(this->button(), &QObject::destroyed, this, [this](){
        if (pressedTimer) pressedTimer->stop();
    });
}

ButtonElement::ButtonElement(const std::string& label, std::function<void ()> onClick)
    : UIElement(new QPushButton(QString::fromStdString(label)))
{
    this->pressedTimer = new QTimer(this);
    this->setOnClick(onClick);
    QObject::connect(this->button(), &QObject::destroyed, this, [this](){
        if (pressedTimer) pressedTimer->stop();
    });
}

QPushButton* ButtonElement::button() {
    return qobject_cast<QPushButton*>(getWidget());
}

ButtonElement& ButtonElement::setOnRepeat(std::function<void ()> onRepeatFunction, int delay_ms, int startup_time)
{
    this->pressedTimer->setInterval(delay_ms);
    this->pressedTimer->setSingleShot(true);
    QObject::connect(this->pressedTimer, &QTimer::timeout, this, [=](){
        if(!(this->button() && this->button()->isDown())) return;
        onRepeatFunction();
        this->pressedTimer->start(delay_ms);
    });
    this->setOnPressed([=]() {
        this->pressedTimer->start(startup_time);
    });
    this->setOnRelease([=]() {
        this->pressedTimer->stop();
    });
    return *this;
}


SliderElement::SliderElement(const std::string& label, float valMin, float valMax, float multiplier, UIElement::LAYOUT orientation)
    : InterfaceUI(orientation)
{
    this->_slider = new FancySlider((orientation == HORIZONTAL ? Qt::Horizontal : Qt::Vertical), valMin, valMax, multiplier);
    this->_label = new LabelElement(label);
    this->add({_label, new UIElement(_slider)});

    defaultValue = valMin;

    QObject::connect(this->slider(), &FancySlider::doubleClicked, this, [=]() { this->setValue(defaultValue); });
}

SliderElement::SliderElement(const std::string& label, float valMin, float valMax, float multiplier, float &binded, UIElement::LAYOUT orientation)
    : SliderElement(label, valMin, valMax, multiplier, orientation)
{
    bindTo(binded);
    this->defaultValue = binded;
}

SliderElement::SliderElement(const std::string &label, float valMin, float valMax, float multiplier, std::function<void (float)> onValueChanged, LAYOUT orientation)
    : SliderElement(label, valMin, valMax, multiplier, orientation)
{
    this->setOnValueChanged(onValueChanged);
}

SliderElement::SliderElement(const std::string &label, float valMin, float valMax, float multiplier, float &binded, std::function<void (float)> onValueChanged, LAYOUT orientation)
    : SliderElement(label, valMin, valMax, multiplier, binded, orientation)
{
    this->setOnValueChanged(onValueChanged);
}

FancySlider* SliderElement::slider() const {
    return this->_slider;
}

LabelElement* SliderElement::label() const
{
    return this->_label;
}

SliderElement& SliderElement::bindTo(float &value) {
    _slider->setfValue(value);
    boundVariable = value;
    setOnValueChanged([this](float newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
    return *this;
}

void SliderElement::update()
{
    if (this->boundVariable.has_value())
        slider()->setfValue(*this->boundVariable);
    return UIElement::update();
}



CheckboxElement::CheckboxElement(const std::string& label)
    : UIElement(new QCheckBox(QString::fromStdString(label)))
{
    get()->setProperty("class", "tight-element");
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

CheckboxElement::CheckboxElement(const std::string& label, bool& binded, std::function<void (bool)> onCheck)
    : CheckboxElement(label, binded)
{
    this->setOnChecked(onCheck);
}

QCheckBox* CheckboxElement::checkBox() {
    return qobject_cast<QCheckBox*>(getWidget());
}

CheckboxElement& CheckboxElement::bindTo(bool &value)
{
    checkBox()->setChecked(value);
    boundVariable = value;
    setOnChecked([this](bool newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
    return *this;
}

void CheckboxElement::update()
{
    if (this->boundVariable.has_value())
        checkBox()->setChecked(*this->boundVariable);
    return UIElement::update();
}

InterfaceUI::InterfaceUI(LAYOUT layout, bool tight, const std::string& title)
    : UIElement(new QGroupBox), title(title)
{
    QLayout* myLayout;
    if (layout == VERTICAL) myLayout = new QVBoxLayout;
    else if (layout == HORIZONTAL) myLayout = new QHBoxLayout;
    else if (layout == GRID) myLayout = new QGridLayout;

    if (tight) {
        myLayout->setSpacing(0);
        myLayout->setContentsMargins(0, 0, 0, 0);
    }
    box()->setStyleSheet(".tight-element{ border: none; }");
    box()->setProperty("class", (tight ? "tight-element" : ""));
    if (layout != GRID)
        getWidget()->setLayout(myLayout);
    // this->setTight(tight);
}
InterfaceUI::InterfaceUI(LAYOUT layout, const std::string& title)
    : InterfaceUI(layout, (title.empty()), title)
{ }

InterfaceUI::~InterfaceUI()
{
    for (auto& child : elements) {
        child->deleteLater();
    }
    elements.clear();
}

QGroupBox* InterfaceUI::box() const
{
    return qobject_cast<QGroupBox*>(getWidget());
}


UIElement* InterfaceUI::add(UIElement* element, const std::string& name)
{
    box()->layout()->addWidget(element->get());
    element->setName(name);
    this->elements.push_back(element);
    this->names.push_back(name);
    return element;
}

void InterfaceUI::add(const std::vector<UIElement*>& elements)
{
    for (auto& element : elements)
        this->add(element);
}

InterfaceUI& InterfaceUI::clear()
{
    for (auto& child : this->elements) {// box()->children()) {
        if (auto interfaceChild = dynamic_cast<InterfaceUI*>(child))
            interfaceChild->clear();
        child->cleanupConnections();
        child->deleteLater();
    }
    this->elements.clear();
    this->names.clear();
    return *this;
}

UIElement* InterfaceUI::findByName(std::string name, bool recursive) const
{
    name = toLower(name);
    for (auto& child : elements) {
        if (toLower(child->getName()) == name) return child;
        if (recursive) {
            if (auto child_as_interface = dynamic_cast<InterfaceUI*>(child)) {
                auto found = child_as_interface->findByName(name, true);
                if (found) return found;
            }
        }
    }
    return nullptr;
}

InterfaceUI *InterfaceUI::setTight(bool tight)
{
    get()->setProperty("class", (tight ? "tight-element" : ""));
}

void InterfaceUI::update()
{
    for (auto& child : this->elements) {// box()->children()) {
        child->update();
    }
    return UIElement::update();
}

InterfaceUI* createHorizontalGroup(std::vector<UIElement*> widgets, bool tight)
{
    auto interface = new InterfaceUI(InterfaceUI::HORIZONTAL, tight);
    for (auto& w : widgets) {
        interface->add(std::move(w));
    }
    return interface;
}

InterfaceUI* createVerticalGroup(std::vector<UIElement*> widgets, bool tight)
{
    auto interface = new InterfaceUI(InterfaceUI::VERTICAL, tight);
    for (auto& w : widgets) {
        interface->add(std::move(w));
    }
    return interface;
}

InterfaceUI* createMultiColumnGroup(std::vector<UIElement*> widgets, int nbColumns, bool tight)
{
    QGridLayout* layout = new QGridLayout;
    for (size_t i = 0; i < widgets.size(); i++) {
        auto& w = widgets[i];
        int row = int (i / nbColumns);
        int col = i % nbColumns;

        layout->addWidget(w->get(), row, col);
    }
    auto ui = new InterfaceUI(InterfaceUI::GRID, tight);
    ui->getWidget()->setLayout(layout);
    return ui;
}

TextEditElement::TextEditElement(const std::string& text, std::string label)
    : InterfaceUI(HORIZONTAL) //UIElement(new QGroupBox)
{
    this->_lineEdit = new QLineEdit(QString::fromStdString(text));
    this->_label = new LabelElement(label);

    /*QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty())
        layout->addWidget(_label);
    layout->addWidget(_lineEdit);
    getWidget()->setLayout(layout);*/
    this->add({_label, new UIElement(_lineEdit)});
}

TextEditElement::TextEditElement(const std::string& text, std::string label, std::string& binded)
    : TextEditElement(text, label)
{
    this->bindTo(binded);
}

QLineEdit* TextEditElement::lineEdit()
{
    return dynamic_cast<QLineEdit*>(_lineEdit);
}

TextEditElement& TextEditElement::setText(const std::string &newText)
{
    lineEdit()->setText(QString::fromStdString(newText));
    return *this;
}

TextEditElement& TextEditElement::setOnTextChange(std::function<void (const std::string&)> func)
{
//    this->addConnection()
    QObject::connect(_lineEdit, &QLineEdit::textChanged, this, [=](QString newText){ // /!\ Capture function by value
        func(newText.toStdString());
    });
    return *this;
}

TextEditElement& TextEditElement::bindTo(std::string& value)
{
    lineEdit()->setText(QString::fromStdString(value));
    boundVariable = value;
    setOnTextChange([this](const std::string& newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
    return *this;
}

void TextEditElement::update()
{
    lineEdit()->setText(QString::fromStdString(*this->boundVariable));
    return UIElement::update();
}


AngleElement::AngleElement(const std::string& label)
    : InterfaceUI(HORIZONTAL) // UIElement(new QGroupBox)
{
    this->_dial = new QDial();
    this->_dial->setWrapping(true); this->_dial->setMinimum(0); this->_dial->setMaximum(360);
    this->_label = new LabelElement(label);

    /*QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty())
        layout->addWidget(_label);
    layout->addWidget(_dial);
    getWidget()->setLayout(layout);*/
    this->add({_label, new UIElement(_dial)});
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

AngleElement& AngleElement::setAngle(float newAngle)
{
    dial()->setValue((int) newAngle);
    if (boundVariable.has_value())
        boundVariable->get() = newAngle;
    return *this;
}

AngleElement& AngleElement::bindTo(float &value)
{
    dial()->setValue((int)value);
    boundVariable = value;
    setOnValueChanged([this](float newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
    return *this;
}

void AngleElement::update()
{
    this->dial()->setValue(*boundVariable);
    return UIElement::update();
}


RangeSliderElement::RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, UIElement::LAYOUT orientation)
    : InterfaceUI(orientation) //UIElement(new QGroupBox)
{
    this->_slider = new RangeSlider((orientation == HORIZONTAL ? Qt::Horizontal : Qt::Vertical), valMin, valMax, multiplier);
    this->_label = new LabelElement(label);

    /*QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty()) {
        layout->addWidget(_label);
    }
    layout->addWidget(_slider);
    getWidget()->setLayout(layout);*/
    this->add({_label, new UIElement(_slider)});
}

RangeSliderElement::RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, float &bindedMin, float &bindedMax, UIElement::LAYOUT orientation)
    : RangeSliderElement(label, valMin, valMax, multiplier, orientation)
{
    this->bindTo(bindedMin, bindedMax);
}

RangeSliderElement::RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, std::function<void(float, float)> onChange, UIElement::LAYOUT orientation)
    : RangeSliderElement(label, valMin, valMax, multiplier, orientation)
{
    this->setOnValueChanged(onChange);
}
RangeSliderElement::RangeSliderElement(const std::string& label, float valMin, float valMax, float multiplier, float& bindedMin, float& bindedMax, std::function<void(float, float)> onChange, UIElement::LAYOUT orientation)
    : RangeSliderElement(label, valMin, valMax, multiplier, bindedMin, bindedMax, orientation)
{
    this->setOnValueChanged(onChange);
}

RangeSlider* RangeSliderElement::slider()
{
    return _slider;
}

LabelElement *RangeSliderElement::label()
{
    return _label;
}

RangeSliderElement& RangeSliderElement::bindTo(float &valueMin, float &valueMax)
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
    return *this;
}

void RangeSliderElement::update()
{
    slider()->setMinValue(*this->boundVariableMin);
    slider()->setMaxValue(*this->boundVariableMax);
    return UIElement::update();
}



ColorPickerElement::ColorPickerElement(const std::string& label)
    : InterfaceUI() // UIElement(new QGroupBox)
{
    this->_colorPicker = new QtColorPicker(this->get());
    this->_label = new LabelElement(label);

    /*
    QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty())
        layout->addWidget(_label);
    layout->addWidget(_colorPicker);
    getWidget()->setLayout(layout);
    */
    this->add({_label, new UIElement(_colorPicker)});
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

ColorPickerElement& ColorPickerElement::setOnSelectionChanged(std::function<void (const Vector3&)> func)
{
    QObject::connect(_colorPicker, &QtColorPicker::colorChanged, this, [=](QColor color) {
        func(ColorPickerElement::qColorToVec3(color));
    });
    return *this;
}

ColorPickerElement& ColorPickerElement::bindTo(Vector3& colorSelected)
{
    boundColor = colorSelected;
    colorPicker()->setCurrentColor(ColorPickerElement::vec3ToQColor(colorSelected));
    this->setOnSelectionChanged([=](Vector3 color) {
        this->boundColor->get() = color;
    });
    return *this;
}

Vector3 ColorPickerElement::getSelection() const
{
    return ColorPickerElement::qColorToVec3(this->colorPicker()->currentColor());
}

void ColorPickerElement::update()
{
    colorPicker()->setCurrentColor(ColorPickerElement::vec3ToQColor(*this->boundColor));
    return UIElement::update();
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

HierarchicalListUI& HierarchicalListUI::setSelectionMode(QAbstractItemView::SelectionMode mode)
{
    this->hierarchicalList()->setSelectionMode(mode);
    return *this;
}


HierarchicalListUI& HierarchicalListUI::clear()
{
    this->hierarchicalList()->clear();
    return *this;
}

HierarchicalListUI& HierarchicalListUI::addItem(HierarchicalListWidgetItemBase* item)
{
    this->hierarchicalList()->addItem(item);
    return *this;
}

HierarchicalListUI& HierarchicalListUI::setCurrentItems(std::vector<int> selectedIDs)
{
    this->hierarchicalList()->setCurrentItems(selectedIDs);
    return *this;
}

std::vector<HierarchicalListWidgetItemBase*> HierarchicalListUI::selectedItems()
{
    auto qItems = this->hierarchicalList()->selectedItems();
    std::vector<HierarchicalListWidgetItemBase*> items;

    for (int i = 0; i < qItems.size(); i++) {
        items.push_back(dynamic_cast<HierarchicalListWidgetItemBase*>(qItems[i]));
    }
    return items;
}

FloatInputElement::FloatInputElement(const std::string &label)
    : InterfaceUI(HORIZONTAL) //UIElement(new QGroupBox)
{
    this->_spinbox = new QDoubleSpinBox(this->get());
    this->_label = new LabelElement(label);

    /*QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty())
        layout->addWidget(_label);
    layout->addWidget(_spinbox);
    getWidget()->setLayout(layout);
    */
    _spinbox->setDecimals(5);
    _spinbox->setMaximum(100000);
    _spinbox->setMinimum(-100000);

    this->add({_label, new UIElement(_spinbox)});
}

FloatInputElement::FloatInputElement(const std::string &label, float &binded)
    : FloatInputElement(label)
{
    this->bindTo(binded);
}

FloatInputElement::FloatInputElement(const std::string& label, std::function<void (float)> onChange)
    : FloatInputElement(label)
{
    QObject::connect(_spinbox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [=](double newValue){ onChange(newValue); });
}

FloatInputElement& FloatInputElement::setValue(float newValue)
{
    spinbox()->setValue(newValue);
    return *this;
}

FloatInputElement& FloatInputElement::bindTo(float& value)
{
    boundVariable = value;
    spinbox()->setValue(value);
    this->setOnValueChanged([=](float newValue) {
        this->boundVariable->get() = newValue;
    });
    return *this;
}

void FloatInputElement::update()
{
    if (this->boundVariable)
        this->spinbox()->setValue(*boundVariable);
    return UIElement::update();
}
