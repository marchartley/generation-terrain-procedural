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



InterfaceUI *PainterToolsUI::createKelvinletToolsUI(ChartView *chartView, PlotModel *dataModel, KelvinletToolParams *params)
{
    params->temporaryVectorData = dataModel->vectorData;

    InterfaceUI* UI = new InterfaceUI(new QVBoxLayout(), true, "Painting tools");

    auto radialScaleSlider = new SliderElement("Radial Scale", params->minRadialScale, params->maxRadialScale, 1.f, params->radialScale);
    auto muSlider = new SliderElement("Mu", params->minMu, params->maxMu, 0.01f, params->mu);
    auto poissonSlider = new SliderElement("Poisson ratio", params->minPoisson, params->maxPoisson, 0.01f, params->poisson);
    auto scaleSlider = new SliderElement("Scale", params->minScale, params->maxScale, 0.01f, params->scale);

    if (!params->currentKelvinlet)
        params->currentKelvinlet = new PinchKelvinlet();
    auto pinchRadio = (new RadioButtonElement("Pinch", [=](bool checked) { if (params->currentKelvinlet) delete params->currentKelvinlet; params->currentKelvinlet = new PinchKelvinlet(); }))->setChecked(typeid(*(params->currentKelvinlet)) == typeid(PinchKelvinlet));
    auto twistRadio = (new RadioButtonElement("Twist", [=](bool checked) { if (params->currentKelvinlet) delete params->currentKelvinlet; params->currentKelvinlet = new TwistKelvinlet(); }))->setChecked(typeid(*(params->currentKelvinlet)) == typeid(TwistKelvinlet));
    auto translateRadio = (new RadioButtonElement("Grab", [=](bool checked) { if (params->currentKelvinlet) delete params->currentKelvinlet; params->currentKelvinlet = new GrabKelvinlet(); }))->setChecked(typeid(*(params->currentKelvinlet)) == typeid(GrabKelvinlet));
    auto scaleRadio = (new RadioButtonElement("Scale", [=](bool checked) { if (params->currentKelvinlet) delete params->currentKelvinlet; params->currentKelvinlet = new ScaleKelvinlet(); }))->setChecked(typeid(*(params->currentKelvinlet)) == typeid(ScaleKelvinlet));

    auto kelvinletsHistory = new HierarchicalListUI();
    // kelvinletsHistory->hierarchicalList()->setSelectionMode();
    for (size_t i = 0; i < params->kelvinlets.size(); i++) {
        const auto& k = params->kelvinlets[i];
        kelvinletsHistory->addItem(new HierarchicalListWidgetItem<Kelvinlet*>(k->getShortName(), k));
    }

    auto editKelvinletCheck = (new CheckboxElement("Edit this Kelvinlet?"))->setChecked(false);
    auto deleteKelvinletButton = new ButtonElement("Delete", [=]() {
        auto items = kelvinletsHistory->selectedItems();
        for (auto& _selection : items) {
            auto selection = dynamic_cast<HierarchicalListWidgetItem<Kelvinlet*>*>(_selection);
            auto found = std::find(params->kelvinlets.begin(), params->kelvinlets.end(), selection->itemValue);
            if (found != params->kelvinlets.end()) {
                params->kelvinlets.erase(found);
            }
            kelvinletsHistory->removeItem(selection->itemValue);
            kelvinletsHistory->update();
        }
    });
    auto validationButton = new ButtonElement("Validate", [=]() {
        GridV3& resultingVectorField = dataModel->vectorData.field;
        resultingVectorField.iterateParallel([&](const Vector3i& p) {
            for (const auto& k : params->kelvinlets)
                resultingVectorField[p] += k->evaluate(p);
        });
        params->kelvinlets.clear();
        kelvinletsHistory->clear();
    });

    UI->add(std::vector<UIElement*>{
        radialScaleSlider,
        muSlider,
        poissonSlider,
        scaleSlider,
        createHorizontalGroupUI(std::vector<UIElement*>{pinchRadio, twistRadio, translateRadio, scaleRadio}),
        kelvinletsHistory,
        createHorizontalGroupUI(std::vector<UIElement*>{editKelvinletCheck, deleteKelvinletButton}),
        validationButton
    });

    QObject::connect(chartView, &ChartView::clickedOnValue, chartView, [=](const Vector3& pos, bool leftClick, bool rightClick) {
        if (pos.isValid()) {
            if (rightClick) {
                params->kelvinletPosition = Vector3::invalid();
            } else {
                params->kelvinletPosition = pos * dataModel->vectorData.getField().getDimensions();
            }
            params->currentKelvinlet->pos = params->kelvinletPosition;
            params->currentKelvinlet->mu = params->mu;
            params->currentKelvinlet->v = params->poisson;
            params->currentKelvinlet->radialScale = params->radialScale;
        }

        if (!pos.isValid() && params->currentKelvinlet->pos.isValid()) {
            const auto& k = params->currentKelvinlet;
            auto clone = k->clone();
            params->kelvinlets.push_back(clone);
            kelvinletsHistory->addItem(new HierarchicalListWidgetItem<Kelvinlet*>(clone->getShortName(), clone));
        }
    });
    QObject::connect(chartView, &ChartView::mouseMoved, chartView, [=](const Vector3& _pos, const Vector3& prevPos) {
        auto pos = _pos * dataModel->vectorData.getField().getDimensions();
        Vector3 force = pos - params->currentKelvinlet->pos;
        if (auto asPinch = dynamic_cast<PinchKelvinlet*>(params->currentKelvinlet)) {
            asPinch->force = force;
        } else if (auto asTwist = dynamic_cast<TwistKelvinlet*>(params->currentKelvinlet)) {
            asTwist->force = Vector3(0, 0, sign(force.x()) * force.norm());
        } else if (auto asGrab = dynamic_cast<GrabKelvinlet*>(params->currentKelvinlet)) {
            asGrab->force = force;
        } else if (auto asScale = dynamic_cast<ScaleKelvinlet*>(params->currentKelvinlet)) {
            asScale->force = sign(force.x()) * force.norm();
        }
    });

    return UI;
}

GridV3& PainterToolsUI::paintKelvinlet(GridV3 &src, const Vector3 &pos, KelvinletToolParams *params)
{
    if (params->currentKelvinlet && params->currentKelvinlet->pos.isValid()) {
        src.iterateParallel([&](const Vector3& p) {
            src[p] += params->currentKelvinlet->evaluate(p);
        });
    }
    return src;
}
