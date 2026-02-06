#include "PainterToolsUI.h"

#include "Utils/Utils.h"


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
        auto grayPicker = new SliderElement("Color", params->minColor.x(), params->maxColor.x(), 0.01f, params->color.x());
        UI->add(grayPicker);
    }

    UI->add(std::vector<UIElement*>{
        createHorizontalGroupUI({addCheck, replaceCheck})
    });

    return UI;
}

GridF getBrushShape(PainterToolParams params) {
    GridF mask(params.radius * 2 + 1, params.radius * 2 + 1, 1);
    if (params.falloff == 0) {
        mask.iterateParallel([&](const Vector3& p) {
            mask[p] = (p - mask.getDimensions().xy() / 2.f).norm2() >= params.radius * params.radius ? 0.f : 1.f;
        });
    } else {
        mask.iterateParallel([&](const Vector3& p) {
            float d = (p - mask.getDimensions().xy() / 2.f).norm();
            float x = d - params.radius * (1.f - params.falloff);
            x /= params.radius * params.falloff;
            mask[p] = 1.f - interpolation::smooth(x, 0.f, 1.f);
        });
    }
    return mask;
}

GridV3& PainterToolsUI::paintImage(GridV3 &src, const Vector3 &pos, PainterToolParams params, bool removeAmount)
{
    GridF mask = getBrushShape(params) * params.opacity;
    // GridV3 maskV3(mask.getDimensions()); maskV3.iterateParallel([&](size_t i){ maskV3[i] = mask[i] * params.color; });
    // src.add(maskV3 * (removeAmount ? -1.f : 1.f), pos - mask.getDimensions().xy() / 2);
    // src.iterateParallel([&](size_t i) {
        // src[i].clamp(params.minClampColor, params.maxClampColor);
    // });
    if (params.additiveMode) {
        GridV3 maskV3(mask.getDimensions()); maskV3.iterateParallel([&](size_t i){ maskV3[i] = mask[i] * params.color; });
        src.add(maskV3 * (removeAmount ? -1.f : 1.f), pos - mask.getDimensions().xy() / 2);
        src.iterateParallel([&](size_t i) {
            src[i].clamp(params.minClampColor, params.maxClampColor);
        });
    } else {
        Vector3i topLeft = pos - mask.getDimensions().xy() / 2;
        mask.iterateParallel([&](const Vector3i& p){
            const float maskVal = mask[p];
            src[topLeft + p] = src[topLeft + p] * (1.f - maskVal) + params.color * maskVal;
            src[topLeft + p].clamp(params.minClampColor, params.maxClampColor);
        });
    }
    return src;
}

GridF& PainterToolsUI::paintImage(GridF &src, const Vector3 &pos, PainterToolParams params, bool removeAmount)
{
    GridF mask = getBrushShape(params) * params.opacity; //GridF::normalizedGaussian(2 * params.radius + 1, 2 * params.radius + 1, 1, float(params.radius * (1.f - params.falloff))) * params.opacity;
    if (params.additiveMode) {
        src.add(mask * params.color.r() * (removeAmount ? -1.f : 1.f), pos - mask.getDimensions().xy() / 2);
        src.iterateParallel([&](size_t i) {
            src[i] = clamp(src[i], params.minClampColor.r(), params.maxClampColor.r());
        });
    } else {
        Vector3i topLeft = pos - mask.getDimensions().xy() / 2;
        mask.iterateParallel([&](const Vector3i& p){
            const float maskVal = mask[p];
            src[topLeft + p] = src[topLeft + p] * (1.f - maskVal) + params.color.r() * maskVal;
            src[topLeft + p] = clamp(src[topLeft + p], params.minClampColor.x(), params.maxClampColor.x());
        });
    }
    return src;
}

