#include "ImageViewer.h"

ImageViewer::ImageViewer(std::string name, QWidget *parent) : ImageViewer(name, new ChartView(new Chart()), parent)
{
}

ImageViewer::ImageViewer(std::string name, ChartView *chartView, QWidget *parent) : AbstractPlotter(name, chartView, parent)
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

    CheckboxElement* displayRButton = (new CheckboxElement(""))->setChecked(this->dataModel->imageData.displayedColors.x());
    CheckboxElement* displayGButton = (new CheckboxElement(""))->setChecked(this->dataModel->imageData.displayedColors.y());
    CheckboxElement* displayBButton = (new CheckboxElement(""))->setChecked(this->dataModel->imageData.displayedColors.z());

    displayRButton->setOnChecked([&](bool toggled) { this->dataModel->imageData.displayedColors.x() = (toggled ? 1 : 0); this->draw(); });
    displayGButton->setOnChecked([&](bool toggled) { this->dataModel->imageData.displayedColors.y() = (toggled ? 1 : 0); this->draw(); });
    displayBButton->setOnChecked([&](bool toggled) { this->dataModel->imageData.displayedColors.z() = (toggled ? 1 : 0); this->draw(); });

    RangeSliderElement* rangeR = (new RangeSliderElement("", 0, 1, 0.01f, Qt::Vertical))->setMinMaxValues(0, 1);
    RangeSliderElement* rangeG = (new RangeSliderElement("", 0, 1, 0.01f, Qt::Vertical))->setMinMaxValues(0, 1);
    RangeSliderElement* rangeB = (new RangeSliderElement("", 0, 1, 0.01f, Qt::Vertical))->setMinMaxValues(0, 1);

    rangeR->setOnValueChanged([&](float newMin, float newMax) { this->dataModel->imageData.colorRangeMin.x() = newMin; this->dataModel->imageData.colorRangeMax.x() = newMax; this->draw(); });
    rangeG->setOnValueChanged([&](float newMin, float newMax) { this->dataModel->imageData.colorRangeMin.y() = newMin; this->dataModel->imageData.colorRangeMax.y() = newMax; this->draw(); });
    rangeB->setOnValueChanged([&](float newMin, float newMax) { this->dataModel->imageData.colorRangeMin.z() = newMin; this->dataModel->imageData.colorRangeMax.z() = newMax; this->draw(); });

    this->viewOptionsInterface->add(std::vector<UIElement*>{
        normalizeModeButton,
        absoluteModeButton,
        createHorizontalGroupUI(std::vector<UIElement*>{
            createVerticalGroupUI(std::vector<UIElement*>({new LabelElement("R"), displayRButton, rangeR})),
            createVerticalGroupUI(std::vector<UIElement*>({new LabelElement("G"), displayGButton, rangeG})),
            createVerticalGroupUI(std::vector<UIElement*>({new LabelElement("B"), displayBButton, rangeB}))
        }),
        createHorizontalGroupUI(std::vector<UIElement*>({
            displayRButton,
            displayGButton,
            displayBButton}))
    });
}

ImageViewer *ImageViewer::getInstance(std::string name)
{
    if (name == "") name = ImageViewer::defaultName;
    if (ImageViewer::instances.count(name) == 0) {
        //        std::cerr << "ImageViewer has not been initialized with function ImageViewer::init()" << std::endl;
        ImageViewer::instances[name] = ImageViewer::init(name);
    }
    return dynamic_cast<ImageViewer*>(ImageViewer::instances[name]);
}

ImageViewer *ImageViewer::init(std::string name, ChartView *chartView, QWidget *parent)
{
    if (ImageViewer::instances.count(name))
        delete ImageViewer::instances[name];
    ImageViewer::instances[name] = new ImageViewer(name, chartView, parent);
    return ImageViewer::getInstance(name);
}

ImageViewer* ImageViewer::updateUI()
{
    AbstractPlotter::updateUI();
    return this;
}

ImageViewer* ImageViewer::setNormalizedModeImage(bool normalize)
{
    this->dataModel->imageData.setNormalized(normalize);
    return updateUI();
}

ImageViewer* ImageViewer::setAbsoluteModeImage(bool absolute)
{
    this->dataModel->imageData.setAbsolute(absolute);
    return updateUI();
}

ImageViewer *ImageViewer::setFilteredValuesImage(bool filtered)
{
    this->dataModel->imageData.setClamped(filtered);
    return updateUI();
}
