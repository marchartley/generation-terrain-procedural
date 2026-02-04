#include "ImageViewerOptionsUI.h"

// ImageViewerOptionsUI::ImageViewerOptionsUI() {}

InterfaceUI *ImageViewerOptionsUI::createRGBImageViewerOptions(ChartView* chartView, PlotModel* dataModel)
{
    InterfaceUI* UI = new InterfaceUI(new QVBoxLayout(), true, "View options");

    auto normalizeModeButton = (new CheckboxElement("Normalize"))->setChecked(false);
    auto absoluteModeButton = (new CheckboxElement("Absolute"))->setChecked(false);

    absoluteModeButton->setOnChecked([=](bool toggled) {
        dataModel->imageData.setAbsolute(toggled);
        chartView->setPlotModel(dataModel);
    });
    normalizeModeButton->setOnChecked([=](bool toggled) {
        dataModel->imageData.setNormalized(toggled);
        chartView->setPlotModel(dataModel);
    });

    CheckboxElement* displayRButton = (new CheckboxElement(""))->setChecked(dataModel->imageData.displayParameters.displayedColors.x());
    CheckboxElement* displayGButton = (new CheckboxElement(""))->setChecked(dataModel->imageData.displayParameters.displayedColors.y());
    CheckboxElement* displayBButton = (new CheckboxElement(""))->setChecked(dataModel->imageData.displayParameters.displayedColors.z());

    displayRButton->setOnChecked([=](bool toggled) { dataModel->imageData.displayParameters.displayedColors.x() = (toggled ? 1 : 0); chartView->setPlotModel(dataModel); });
    displayGButton->setOnChecked([=](bool toggled) { dataModel->imageData.displayParameters.displayedColors.y() = (toggled ? 1 : 0); chartView->setPlotModel(dataModel); });
    displayBButton->setOnChecked([=](bool toggled) { dataModel->imageData.displayParameters.displayedColors.z() = (toggled ? 1 : 0); chartView->setPlotModel(dataModel); });

    RangeSliderElement* rangeR = (new RangeSliderElement("", 0, 1, 0.01f, Qt::Vertical))->setMinMaxValues(0, 1);
    RangeSliderElement* rangeG = (new RangeSliderElement("", 0, 1, 0.01f, Qt::Vertical))->setMinMaxValues(0, 1);
    RangeSliderElement* rangeB = (new RangeSliderElement("", 0, 1, 0.01f, Qt::Vertical))->setMinMaxValues(0, 1);

    rangeR->setOnValueChanged([=](float newMin, float newMax) { dataModel->imageData.displayParameters.colorRangeMin.x() = newMin; dataModel->imageData.displayParameters.colorRangeMax.x() = newMax; chartView->setPlotModel(dataModel); });
    rangeG->setOnValueChanged([=](float newMin, float newMax) { dataModel->imageData.displayParameters.colorRangeMin.y() = newMin; dataModel->imageData.displayParameters.colorRangeMax.y() = newMax; chartView->setPlotModel(dataModel); });
    rangeB->setOnValueChanged([=](float newMin, float newMax) { dataModel->imageData.displayParameters.colorRangeMin.z() = newMin; dataModel->imageData.displayParameters.colorRangeMax.z() = newMax; chartView->setPlotModel(dataModel); });

    std::vector<UIElement*> overlayCheckboxes;
    for (auto& [name, over] : chartView->overlayColors) {
        CheckboxElement* checkOverlay = new CheckboxElement(name, [=](bool checked) { chartView->overlayDisplayed[name] = checked; chartView->setPlotModel(dataModel); });
        checkOverlay->setChecked(chartView->overlayDisplayed[name]);
        overlayCheckboxes.push_back(checkOverlay);
    }

    UI->add(std::vector<UIElement*>{
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
                                                          displayBButton})),
        createVerticalGroupUI(overlayCheckboxes)
    });
    return UI;
}

InterfaceUI *ImageViewerOptionsUI::createGreyImageViewerOptions(ChartView *chartView, PlotModel *dataModel)
{
    InterfaceUI* UI = new InterfaceUI(new QVBoxLayout(), true, "View options");

    auto normalizeModeButton = (new CheckboxElement("Normalize"))->setChecked(false);
    auto absoluteModeButton = (new CheckboxElement("Absolute"))->setChecked(false);

    absoluteModeButton->setOnChecked([=](bool toggled) {
        dataModel->imageData.setAbsolute(toggled);
        chartView->setPlotModel(dataModel);
    });
    normalizeModeButton->setOnChecked([=](bool toggled) {
        dataModel->imageData.setNormalized(toggled);
        chartView->setPlotModel(dataModel);
    });

    float mini = dataModel->getImage().min().r();
    float maxi = dataModel->getImage().max().r();
    RangeSliderElement* rangeSlider = (new RangeSliderElement("", mini, maxi, 0.01f, Qt::Vertical))->setMinMaxValues(mini, maxi);

    rangeSlider->setOnValueChanged([=](float newMin, float newMax) {
        dataModel->imageData.displayParameters.colorRangeMin.x() = newMin;
        // dataModel->imageData.displayParameters.colorRangeMin.y() = newMin;
        // dataModel->imageData.displayParameters.colorRangeMin.z() = newMin;
        dataModel->imageData.displayParameters.colorRangeMax.x() = newMax;
        // dataModel->imageData.displayParameters.colorRangeMax.y() = newMax;
        // dataModel->imageData.displayParameters.colorRangeMax.z() = newMax;
        chartView->setPlotModel(dataModel);
    });

    std::vector<UIElement*> overlayCheckboxes;
    for (auto& [name, over] : chartView->overlayColors) {
        CheckboxElement* checkOverlay = new CheckboxElement(name, [=](bool checked) { chartView->overlayDisplayed[name] = checked; chartView->setPlotModel(dataModel); });
        checkOverlay->setChecked(chartView->overlayDisplayed[name]);
        overlayCheckboxes.push_back(checkOverlay);
    }

    UI->add(std::vector<UIElement*>{
        normalizeModeButton,
        absoluteModeButton,
        createVerticalGroupUI(std::vector<UIElement*>({new LabelElement("Value"), rangeSlider})),
        createVerticalGroupUI(overlayCheckboxes)
    });
    return UI;
}
