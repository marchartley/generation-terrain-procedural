#include "ImageViewerOptionsUI.h"

// ImageViewerOptionsUI::ImageViewerOptionsUI() {}

InterfaceUI* ImageViewerOptionsUI::createRGBImageViewerOptions(AbstractPlotter* plotter)
{
    // auto& dataModel = plotter->dataModel;
    // auto& chartView = plotter->chartView;
    const GridV3 img = plotter->dataModel->imageData.image.getColorImage();

    auto UI = new InterfaceUI(new QVBoxLayout(), true, "View options");

    auto normalizeModeButton = new CheckboxElement("Normalize");
    auto absoluteModeButton = new CheckboxElement("Absolute");

    auto displayRButton = new CheckboxElement("");
    auto displayGButton = new CheckboxElement("");
    auto displayBButton = new CheckboxElement("");

    auto rangeR = new RangeSliderElement("", 0, 1, 0.01f);
    auto rangeG = new RangeSliderElement("", 0, 1, 0.01f);
    auto rangeB = new RangeSliderElement("", 0, 1, 0.01f);

    normalizeModeButton->setChecked(false);
    absoluteModeButton->setChecked(false);

    absoluteModeButton->setOnChecked([=](bool toggled) {
        plotter->dataModel->imageData.setAbsolute(toggled);
        plotter->chartView->setPlotModel(plotter->dataModel);
    });
    normalizeModeButton->setOnChecked([=](bool toggled) {
        plotter->dataModel->imageData.setNormalized(toggled);
        plotter->chartView->setPlotModel(plotter->dataModel);
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


    displayRButton->setChecked(plotter->dataModel->imageData.displayParameters.displayedColors.x()).setOnChecked([=](bool toggled) { plotter->dataModel->imageData.displayParameters.displayedColors.x() = (toggled ? 1 : 0); plotter->chartView->setPlotModel(plotter->dataModel); });
    displayGButton->setChecked(plotter->dataModel->imageData.displayParameters.displayedColors.y()).setOnChecked([=](bool toggled) { plotter->dataModel->imageData.displayParameters.displayedColors.y() = (toggled ? 1 : 0); plotter->chartView->setPlotModel(plotter->dataModel); });
    displayBButton->setChecked(plotter->dataModel->imageData.displayParameters.displayedColors.z()).setOnChecked([=](bool toggled) { plotter->dataModel->imageData.displayParameters.displayedColors.z() = (toggled ? 1 : 0); plotter->chartView->setPlotModel(plotter->dataModel); });


    rangeR->setMinMaxValues(0, 1).setOnValueChanged([=](float newMin, float newMax) { plotter->dataModel->imageData.displayParameters.colorRangeMin.x() = newMin; plotter->dataModel->imageData.displayParameters.colorRangeMax.x() = newMax; plotter->chartView->setPlotModel(plotter->dataModel); });
    rangeG->setMinMaxValues(0, 1).setOnValueChanged([=](float newMin, float newMax) { plotter->dataModel->imageData.displayParameters.colorRangeMin.y() = newMin; plotter->dataModel->imageData.displayParameters.colorRangeMax.y() = newMax; plotter->chartView->setPlotModel(plotter->dataModel); });
    rangeB->setMinMaxValues(0, 1).setOnValueChanged([=](float newMin, float newMax) { plotter->dataModel->imageData.displayParameters.colorRangeMin.z() = newMin; plotter->dataModel->imageData.displayParameters.colorRangeMax.z() = newMax; plotter->chartView->setPlotModel(plotter->dataModel); });

    std::vector<UIElement*> overlayCheckboxes;
    for (auto& [name, over] : plotter->chartView->overlayColors) {
        auto checkOverlay = new CheckboxElement(name, [=](bool checked) { plotter->chartView->overlayDisplayed[name] = checked; plotter->chartView->setPlotModel(plotter->dataModel); });
        checkOverlay->setChecked(plotter->chartView->overlayDisplayed[name]);
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

InterfaceUI* ImageViewerOptionsUI::createGreyImageViewerOptions(AbstractPlotter* plotter)
{
    // auto& dataModel = plotter->dataModel;
    // auto& chartView = plotter->chartView;

    auto greyImg = plotter->dataModel->getImageGrey();
    float mini = greyImg.min();
    float maxi = greyImg.max();

    auto UI = new InterfaceUI(new QVBoxLayout(), true, "View options");

    auto normalizeModeButton = new CheckboxElement("Normalize");
    auto absoluteModeButton = new CheckboxElement("Absolute");

    auto rangeSlider = new RangeSliderElement("", mini, maxi, 0.01f, Qt::Vertical);



    normalizeModeButton->setChecked(false);
    absoluteModeButton->setChecked(false);

    absoluteModeButton->setOnChecked([=](bool toggled) {
        plotter->dataModel->imageData.setAbsolute(toggled);
        plotter->chartView->setPlotModel(plotter->dataModel);
    });
    normalizeModeButton->setOnChecked([=](bool toggled) {
        plotter->dataModel->imageData.setNormalized(toggled);
        plotter->chartView->setPlotModel(plotter->dataModel);
    });

    rangeSlider->setMinMaxValues(mini, maxi);

    rangeSlider->setOnValueChanged([=](float newMin, float newMax) {
        plotter->dataModel->imageData.displayParameters.colorRangeMin.x() = newMin;
        plotter->dataModel->imageData.displayParameters.colorRangeMax.x() = newMax;
        plotter->chartView->setPlotModel(plotter->dataModel);
    });

    std::vector<UIElement*> overlayCheckboxes;
    for (auto& [name, over] : plotter->chartView->overlayColors) {
        auto checkOverlay = new CheckboxElement(name, [=](bool checked) { plotter->chartView->overlayDisplayed[name] = checked; plotter->chartView->setPlotModel(plotter->dataModel); });
        checkOverlay->setChecked(plotter->chartView->overlayDisplayed[name]);
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
