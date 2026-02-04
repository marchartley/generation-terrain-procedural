#include "PainterToolsUI.h"



InterfaceUI *PainterToolsUI::createPainterToolsUI(ChartView *chartView, PlotModel *dataModel, PainterToolParams *params)
{
    InterfaceUI* UI = new InterfaceUI(new QVBoxLayout(), true, "Painting tools");

    auto radiusSlider = new SliderElement("Radius", params->minRadius, params->maxRadius, 1.f, params->radius);
    auto amountSlider = new SliderElement("Opacity", params->minOpacity, params->maxOpacity, 0.01f, params->opacity);
    auto falloffSlider = new SliderElement("Blur", 0.f, 1.f, 0.01f, params->falloff);

    auto addCheck = (new RadioButtonElement("Add"))->setOnChecked([=](bool checked){ params->additiveMode = checked; });
    auto replaceCheck = (new RadioButtonElement("Replace"))->setOnChecked([=](bool checked){ params->additiveMode = !checked; });
    addCheck ->setChecked(params->additiveMode);
    replaceCheck->setChecked(!params->additiveMode);

    UI->add(std::vector<UIElement*>{
        radiusSlider,
        amountSlider,
        falloffSlider});
    if (params->RGBimage) {
        auto colorPicker = new ColorPickerElement("Color", params->color);
        UI->add(colorPicker);
    } else {
        auto grayPicker = new SliderElement("Color", 0.f, 1.f, 0.01f, params->color.x());
        UI->add(grayPicker);
    }

    UI->add(std::vector<UIElement*>{
        createHorizontalGroupUI({addCheck, replaceCheck})
    });

    return UI;
}



GridV3& PainterToolsUI::paintImage(GridV3 &src, const Vector3 &pos, PainterToolParams params, bool removeAmount)
{
    GridF mask = GridF::normalizedGaussian(2 * params.radius + 1, 2 * params.radius + 1, 1, float(params.radius * (1.f - params.falloff))) * params.opacity;
    // if (params.RGBimage) {
    GridV3 maskV3(mask.getDimensions()); maskV3.iterateParallel([&](size_t i){ maskV3[i] = mask[i] * params.color; });
    src.add(maskV3 * (removeAmount ? -1.f : 1.f), pos - mask.getDimensions().xy() / 2);
    src.iterateParallel([&](size_t i) {
        src[i].clamp(params.minColor, params.maxColor);
    });
    // } else {
    // src.add(mask * params.color.r() * (removeAmount ? -1.f : 1.f), pos - mask.getDimensions().xy() / 2);
    // }
    // Q_EMIT this->imagePainted(src);
    // this->draw();
}

GridF& PainterToolsUI::paintImage(GridF &src, const Vector3 &pos, PainterToolParams params, bool removeAmount)
{
    GridF mask = GridF::normalizedGaussian(2 * params.radius + 1, 2 * params.radius + 1, 1, float(params.radius * (1.f - params.falloff))) * params.opacity;
    // if (params.RGBimage) {
    // GridF maskF(mask.getDimensions()); maskF.iterateParallel([&](size_t i){ maskF[i] = mask[i] * params.color; });
    // src.add(maskF * (removeAmount ? -1.f : 1.f), pos - mask.getDimensions().xy() / 2);
    // } else {
    src.add(mask * params.color.r() * (removeAmount ? -1.f : 1.f), pos - mask.getDimensions().xy() / 2);
    src.iterateParallel([&](size_t i) {
        src[i] = clamp(src[i], params.minColor.r(), params.maxColor.r());
    });
    // }
    // Q_EMIT this->imagePainted(src);
    // this->draw();
}

