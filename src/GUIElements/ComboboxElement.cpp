#include "ComboboxElement.h"


ComboboxElement::ComboboxElement(const std::string& label)
    : UIElement(new QGroupBox)
{
    this->_combobox = new QComboBox();
    this->_label = new QLabel(QString::fromStdString(label));

    QBoxLayout* layout = new QHBoxLayout;
    layout->setMargin(0);
    if (!label.empty())
        layout->addWidget(_label);
    layout->addWidget(_combobox);
    getWidget()->setLayout(layout);
}

ComboboxElement::ComboboxElement(const std::string& label, std::vector<ComboboxLineElementBase*> choices)
    : ComboboxElement(label)
{
    this->choices = choices;
    for (auto& c : choices) {
        if (c->iconPath != "") {
            combobox()->addItem(QIcon(QString::fromStdString(c->iconPath)), QString::fromStdString(c->label));
            // combobox()->addItem(QIcon(QString::fromStdString(c->iconPath)), QString::fromStdString(c->label), c->label);
        } else {
            combobox()->addItem(QString::fromStdString(c->label));
            // combobox()->addItem(QString::fromStdString(c->label), c->label);
        }
    }
}

ComboboxElement::ComboboxElement(const std::string& label, std::vector<ComboboxLineElementBase*> choices, int &currentSelection)
    : ComboboxElement(label, choices)
{
    this->bindTo(currentSelection);
}

QComboBox* ComboboxElement::combobox() const
{
    return _combobox;
}

ComboboxElement* ComboboxElement::setOnSelectionChanged(std::function<void (int)> func)
{
    QObject::connect(_combobox, &QComboBox::currentTextChanged, this, [=](QString text) {
        int index = _combobox->currentIndex();
        func(index);
    });
    return this;
}

ComboboxElement* ComboboxElement::bindTo(int &indexSelected)
{
    boundIndex = indexSelected;
    combobox()->setCurrentIndex(indexSelected);
    this->setOnSelectionChanged([&](int index) {
        this->boundIndex->get() = index;
    });
    return this;
}

ComboboxLineElementBase* ComboboxElement::getSelection() const
{
    if(combobox()->currentIndex() >= 0)
        return this->choices.at(combobox()->currentIndex());
    return nullptr;
}

ComboboxElement *ComboboxElement::addChoice(ComboboxLineElementBase* addedChoice, bool selected)
{
    this->choices.push_back(addedChoice);
    if (selected && this->boundIndex)
        this->boundIndex->get() = this->choices.size() - 1;
    this->update();
    if (selected)
        combobox()->setCurrentIndex(this->choices.size() - 1);
    return this;
}

void ComboboxElement::update()
{
    int currentSelection = combobox()->currentIndex();
    combobox()->clear();
    for (auto& c : choices) {
        if (c->iconPath != "") {
            combobox()->addItem(QIcon(QString::fromStdString(c->iconPath)), QString::fromStdString(c->label));
            // combobox()->addItem(QIcon(QString::fromStdString(c->iconPath)), QString::fromStdString(c->label), c->label);
        } else {
            combobox()->addItem(QString::fromStdString(c->label));
            // combobox()->addItem(QString::fromStdString(c->label), c->label);
        }
    }
    if (this->boundIndex)
        combobox()->setCurrentIndex(*this->boundIndex);
    else if (currentSelection >= 0 && currentSelection < (int) choices.size())
        combobox()->setCurrentIndex(currentSelection);
}



ComboboxLineElementBase::ComboboxLineElementBase()
    : ComboboxLineElementBase("")
{

}

ComboboxLineElementBase::ComboboxLineElementBase(const std::string& label)
    : ComboboxLineElementBase(label, "")
{

}

ComboboxLineElementBase::ComboboxLineElementBase(const std::string& label, const std::string& iconPath)
    : label(label), iconPath(iconPath)
{

}
