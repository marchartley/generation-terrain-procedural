#include "ImageViewerOptionsUI.h"

// ImageViewerOptionsUI::ImageViewerOptionsUI() {}

InterfaceUI *ImageViewerOptionsUI::createRGBImageViewerOptions(ChartView* chartView, PlotModel* dataModel)
{
    const GridV3 img = dataModel->imageData.image.getColorImage();

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

    Vector3 minColors = Vector3::max();
    Vector3 maxColors = Vector3::min();
    img.iterate([&](size_t i) {
        minColors.r() = std::min(minColors.r(), img[i].r());
        minColors.g() = std::min(minColors.g(), img[i].g());
        minColors.b() = std::min(minColors.b(), img[i].b());
        maxColors.r() = std::max(maxColors.r(), img[i].r());
        maxColors.g() = std::max(maxColors.g(), img[i].g());
        maxColors.b() = std::max(maxColors.b(), img[i].b());
    });

    if (img.empty()) {
        minColors = Vector3(0, 0, 0);
        maxColors = Vector3(0, 0, 0);
    }

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
            createVerticalGroupUI(std::vector<UIElement*>({new LabelElement("R"), new LabelElement(std::to_string(minColors.r())), rangeR, new LabelElement(std::to_string(maxColors.r())), displayRButton})),
            createVerticalGroupUI(std::vector<UIElement*>({new LabelElement("G"), new LabelElement(std::to_string(minColors.g())), rangeG, new LabelElement(std::to_string(maxColors.g())), displayGButton})),
            createVerticalGroupUI(std::vector<UIElement*>({new LabelElement("B"), new LabelElement(std::to_string(minColors.b())), rangeB, new LabelElement(std::to_string(maxColors.b())), displayBButton}))
        }),
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

    auto greyImg = dataModel->getImageGrey();
    float mini = greyImg.min();
    float maxi = greyImg.max();
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
        createHorizontalGroupUI({
            createVerticalGroupUI(std::vector<UIElement*>{new LabelElement("Value"), rangeSlider}),
            createVerticalGroupUI(std::vector<UIElement*>{
                                 new LabelElement("Min: " + std::to_string(mini)),
                                 new LabelElement("Max: " + std::to_string(maxi))
            })
        }),
        createVerticalGroupUI(overlayCheckboxes)
    });
    return UI;
}
