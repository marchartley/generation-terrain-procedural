#include "Plotter.h"

Plotter::Plotter(std::string name, QWidget *parent) : Plotter(name, new ChartView(new Chart()), parent)
{
}

Plotter::Plotter(std::string name, ChartView *chartView, QWidget *parent) : AbstractPlotter(name, chartView, parent)
{
    auto normalizeModeButton = (new CheckboxElement("Normalize"))->setChecked(false);
    auto absoluteModeButton = (new CheckboxElement("Absolute"))->setChecked(false);
    // this->rangeValuesWidget = new RangeSliderElement("Values", -1000, 1000, 0.01f, -1000, 1000, Qt::Vertical);
    this->rangeValuesWidget = (new RangeSliderElement("Values", -1000, 1000, 0.01f, Qt::Vertical))->setMinMaxValues(-1000, 1000);
    // auto rangeActiveCheckbox = (new CheckboxElement("Filter"))->setChecked(false);

    absoluteModeButton->setOnChecked([&](bool toggled) {
        this->dataModel->imageData.setAbsolute(toggled);
        this->draw();
    });
    normalizeModeButton->setOnChecked([&](bool toggled) {
        this->dataModel->imageData.setNormalized(toggled);
        this->draw();
    });
    /*rangeValuesWidget->setOnValueChanged([&](float minVal, float maxVal) {
        this->dataModel->imageData.setColorRanges(minVal, maxVal);
        this->draw();
    });*/

    CheckboxElement* displayRButton = (new CheckboxElement("R"))->setChecked(this->dataModel->imageData.displayedColors.x());
    CheckboxElement* displayGButton = (new CheckboxElement("G"))->setChecked(this->dataModel->imageData.displayedColors.y());
    CheckboxElement* displayBButton = (new CheckboxElement("B"))->setChecked(this->dataModel->imageData.displayedColors.z());

    displayRButton->setOnChecked([&](bool toggled) { this->dataModel->imageData.displayedColors.x() = (toggled ? 1 : 0); this->draw(); });
    displayGButton->setOnChecked([&](bool toggled) { this->dataModel->imageData.displayedColors.y() = (toggled ? 1 : 0); this->draw(); });
    displayBButton->setOnChecked([&](bool toggled) { this->dataModel->imageData.displayedColors.z() = (toggled ? 1 : 0); this->draw(); });

    this->viewOptionsInterface->add(std::vector<UIElement*>({normalizeModeButton, absoluteModeButton, rangeValuesWidget, createHorizontalGroupUI(std::vector<UIElement*>{displayRButton, displayGButton, displayBButton})}));
}

Plotter *Plotter::getInstance(std::string name)
{
    if (name == "") name = Plotter::defaultName;
    if (Plotter::instances.count(name) == 0) {
        //        std::cerr << "Plotter has not been initialized with function Plotter::init()" << std::endl;
        Plotter::instances[name] = Plotter::init(name);
    }
    return dynamic_cast<Plotter*>(Plotter::instances[name]);
}

Plotter *Plotter::init(std::string name, ChartView *chartView, QWidget *parent)
{
    if (Plotter::instances.count(name))
        delete Plotter::instances[name];
    Plotter::instances[name] = new Plotter(name, chartView, parent);
    return Plotter::getInstance(name);
}

Plotter* Plotter::updateUI()
{
    AbstractPlotter::updateUI();
    return this;
}

Plotter* Plotter::setNormalizedModeImage(bool normalize)
{
    this->dataModel->imageData.setNormalized(normalize);
    return updateUI();
}

Plotter* Plotter::setAbsoluteModeImage(bool absolute)
{
    this->dataModel->imageData.setAbsolute(absolute);
    return updateUI();
}

Plotter *Plotter::setFilteredValuesImage(bool filtered)
{
    this->dataModel->imageData.setClamped(filtered);
    return updateUI();
}
