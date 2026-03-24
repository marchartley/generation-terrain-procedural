#ifndef COMBOBOXELEMENT_H
#define COMBOBOXELEMENT_H

#include "GUIElements/CommonInterface.h"
#include <QComboBox>

class ComboboxLineElementBase;
template <class T>
class ComboboxLineElement;

class ComboboxElement : public InterfaceUI {
    Q_OBJECT
public:
    ComboboxElement(const std::string& label);
    ComboboxElement(const std::string& label, std::vector<ComboboxLineElementBase*> choices);
    ComboboxElement(const std::string& label, std::vector<ComboboxLineElementBase*> choices, int& currentSelection);

    template <class T>
    ComboboxElement(const std::string& label, std::vector<ComboboxLineElement<T>*> choices);
    template <class T>
    ComboboxElement(const std::string& label, std::vector<ComboboxLineElement<T>*> choices, int& currentSelection);

    QComboBox* combobox() const;

    ComboboxElement& setOnSelectionChanged(std::function<void(int)> func);

    ComboboxElement& bindTo(int& indexSelected);

    ComboboxLineElementBase* getSelection() const;

    template <class T>
    T getSelection() const;

    ComboboxElement& addChoice(ComboboxLineElementBase* addedChoice, bool selected = false);

public Q_SLOTS:
    void update();

    //protected:
public:
    std::vector<ComboboxLineElementBase*> choices;
    std::optional<std::reference_wrapper<int>> boundIndex;
    //    std::optional<std::reference_wrapper<std::vector<std::string>>> boundValues;

    LabelElement* _label;
    QComboBox* _combobox;

    //    bool itemsAreImages = false;
};


class ComboboxLineElementBase {
public:
    ComboboxLineElementBase();
    ComboboxLineElementBase(const std::string& label);
    ComboboxLineElementBase(const std::string& label, const std::string& iconPath);

    virtual ~ComboboxLineElementBase() {}

    std::string label;
    std::string iconPath;
};

template <class T = int>
class ComboboxLineElement : public ComboboxLineElementBase {
public:
    ComboboxLineElement();
    ComboboxLineElement(const std::string& label);
    ComboboxLineElement(const std::string& label, T value);
    ComboboxLineElement(const std::string& label, const std::string& iconPath, T value);

    T value;
};




template<class T>
T ComboboxElement::getSelection() const
{
    auto res = this->getSelection();
    return dynamic_cast<ComboboxLineElement<T>*>(res)->value;
}

template <class T>
ComboboxElement::ComboboxElement(const std::string &label, std::vector<ComboboxLineElement<T> *> choices)
    : ComboboxElement(label)
{
    auto newChoices = std::vector<ComboboxLineElementBase*>(choices.size());
    for (size_t i = 0; i < choices.size(); i++) {
        auto& choice = choices[i];
        newChoices[i] = dynamic_cast<ComboboxLineElementBase*>(choice);
    }
    this->choices = newChoices;
    for (auto& c : newChoices) {
        if (c->iconPath != "") {
            combobox()->addItem(QIcon(QString::fromStdString(c->iconPath)), QString::fromStdString(c->label));
            // combobox()->addItem(QIcon(QString::fromStdString(c->iconPath)), QString::fromStdString(c->label), c->label);
        } else {
            combobox()->addItem(QString::fromStdString(c->label));
            // combobox()->addItem(QString::fromStdString(c->label), c->label);
        }
    }
}
template <class T>
ComboboxElement::ComboboxElement(const std::string &label, std::vector<ComboboxLineElement<T> *> choices, int& currentSelection)
    : ComboboxElement(label, choices)
{
    this->bindTo(currentSelection);
}



template <class T>
ComboboxLineElement<T>::ComboboxLineElement()
    : ComboboxLineElement("")
{
}

template <class T>
ComboboxLineElement<T>::ComboboxLineElement(const std::string& label)
    : ComboboxLineElement<T>(label, -1)
{

}

template <class T>
ComboboxLineElement<T>::ComboboxLineElement(const std::string& label, T value)
    : ComboboxLineElement<T>(label, "", value)
{

}

template <class T>
ComboboxLineElement<T>::ComboboxLineElement(const std::string &label, const std::string &iconPath, T value)
    : ComboboxLineElementBase(label, iconPath), value(value)
{
}

#endif // COMBOBOXELEMENT_H
