#include "CommonInterface.h"

#include <QBoxLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>

UIElement::UIElement(QWidget *widget)
    : element(widget)
{

}

QWidget *UIElement::getWidget() const {
    return element;
}

void UIElement::cleanupConnections() {
    for (const auto& connection : connections) {
        QObject::disconnect(connection);
    }
    connections.clear();
}


UIElement::~UIElement() {
    cleanupConnections();
    delete element;
}

void UIElement::setName(std::string name)
{
    this->name = name;
}
const std::string &UIElement::getName() const
{
    return name;
}




ButtonElement::ButtonElement(std::string label)
    : UIElement(new QPushButton(QString::fromStdString(label))) {}

ButtonElement::ButtonElement(std::string label, std::function<void ()> onClick)
    : UIElement(new QPushButton(QString::fromStdString(label)))
{
    this->setOnClick(onClick);
}

QPushButton *ButtonElement::button() {
    return static_cast<QPushButton*>(getWidget());
}




SliderElement::SliderElement(std::string label, float valMin, float valMax, float multiplier, Qt::Orientation orientation)
    : UIElement(new QGroupBox)
{
    this->_slider = new FancySlider(orientation, valMin, valMax, multiplier);
    this->_label = new QLabel(QString::fromStdString(label));

    QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    layout->addWidget(_label);
    layout->addWidget(_slider);
    getWidget()->setLayout(layout);
}

SliderElement::SliderElement(std::string label, float valMin, float valMax, float multiplier, float &binded, Qt::Orientation orientation)
    : SliderElement(label, valMin, valMax, multiplier, orientation)
{
    bindTo(binded);
}

FancySlider* SliderElement::slider() {
    return this->_slider;
}

QLabel *SliderElement::label()
{
    return this->_label;
}

void SliderElement::bindTo(float &value) {
    _slider->setfValue(value);
    boundVariable = value;
    setOnValueChanged([this](float newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
}



CheckboxElement::CheckboxElement(std::string label)
    : UIElement(new QCheckBox(QString::fromStdString(label)))
{

}

CheckboxElement::CheckboxElement(std::string label, bool &binded)
    : CheckboxElement(label)
{
    bindTo(binded);
}

CheckboxElement::CheckboxElement(std::string label, std::function<void (bool)> onCheck)
    : CheckboxElement(label)
{
    this->setOnChecked(onCheck);
}

QCheckBox *CheckboxElement::checkBox() {
    return static_cast<QCheckBox*>(getWidget());
}

void CheckboxElement::bindTo(bool &value)
{
    checkBox()->setChecked(value);
    boundVariable = value;
    setOnChecked([this](bool newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
}

InterfaceUI::InterfaceUI(QLayout *layout, std::string title)
    : UIElement(new QGroupBox), title(title)
{
    getWidget()->setLayout(layout);
}

InterfaceUI::~InterfaceUI()
{
    elements.clear();
}

QGroupBox *InterfaceUI::box() const
{
    return static_cast<QGroupBox*>(getWidget());
}

UIElement *InterfaceUI::add(UIElement *element, std::string name)
{
    element->setName(name);
    this->elements.push_back(element);
    this->names.push_back(name);
    box()->layout()->addWidget(element->get());
    return element;
}

UIElement *InterfaceUI::add(QLayout *layout, std::string name)
{
    InterfaceUI* interface = new InterfaceUI(layout, name);
    return this->add(interface, name);
}

UIElement *InterfaceUI::find(std::string name)
{
    for (size_t i = 0; i < names.size(); i++)
        if (names[i] == name)
            return elements[i];
    return nullptr;
}

RadioButtonElement::RadioButtonElement(std::string label)
    : UIElement(new QRadioButton(QString::fromStdString(label)))
{

}

RadioButtonElement::RadioButtonElement(std::string label, bool &binded)
    : RadioButtonElement(label)
{
    bindTo(binded);
}

QRadioButton *RadioButtonElement::radioButton()
{
    return static_cast<QRadioButton*>(getWidget());
}

void RadioButtonElement::bindTo(bool &value)
{
    radioButton()->setChecked(value);
    boundVariable = value;
    setOnChecked([this](bool newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
}

InterfaceUI *createHorizontalGroupUI(std::vector<UIElement *> widgets)
{
    InterfaceUI* interface = new InterfaceUI(new QHBoxLayout);
    for (UIElement*& w : widgets)
        interface->add(w);
    return interface;
}

InterfaceUI *createVerticalGroupUI(std::vector<UIElement *> widgets)
{
    InterfaceUI* interface = new InterfaceUI(new QVBoxLayout);
    for (UIElement*& w : widgets)
        interface->add(w);
    return interface;
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

QLineEdit *TextEditElement::lineEdit()
{
    return dynamic_cast<QLineEdit*>(_lineEdit);
}

void TextEditElement::setOnTextChange(std::function<void (std::string)> func)
{
//    this->addConnection()
    QObject::connect(_lineEdit, &QLineEdit::textChanged, this, [=](QString newText){ // /!\ Capture function by value
        func(newText.toStdString());
    });
}

void TextEditElement::bindTo(std::string &value)
{
    lineEdit()->setText(QString::fromStdString(value));
    boundVariable = value;
    setOnTextChange([this](std::string newValue) {
        if (boundVariable) {
            boundVariable->get() = newValue;
        }
    });
}
