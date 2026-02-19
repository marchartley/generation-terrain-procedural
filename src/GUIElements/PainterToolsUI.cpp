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


void updateKelvinletList(HierarchicalListUI* listWidget, KelvinletToolParams* params) {
    listWidget->get()->blockSignals(true);
    auto currentSelection = listWidget->hierarchicalList()->currentRow();
    listWidget->clear();
    for (size_t i = 0; i < params->kelvinlets.size(); i++) {
        const auto& k = params->kelvinlets[i];
        listWidget->addItem(new HierarchicalListWidgetItem<Kelvinlet*>(k->getShortName(), k));
    }
    for (size_t i = 0; i < params->additional_kelvinlets.size(); i++) {
        const auto& k = params->additional_kelvinlets[i];
        listWidget->addItem(new HierarchicalListWidgetItem<Kelvinlet*>("[*] " + k->getShortName(), k));
    }
    listWidget->hierarchicalList()->setCurrentRow(currentSelection);
    listWidget->get()->blockSignals(false);
    listWidget->update();
}
InterfaceUI *PainterToolsUI::createKelvinletToolsUI(AbstractPlotter* plotter, ChartView* chartView, PlotModel* dataModel, KelvinletToolParams* params)
{
    params->temporaryVectorData = dataModel->vectorData;

    InterfaceUI* UI = new InterfaceUI(new QVBoxLayout(), false, "Kelvinlet tool editor");

    auto radialScaleSlider = new SliderElement("Radial Scale", params->minRadialScale, params->maxRadialScale, 1.f, params->radialScale);
    auto muSlider = new SliderElement("Mu", params->minMu, params->maxMu, 0.01f, params->mu);
    auto poissonSlider = new SliderElement("Poisson ratio", params->minPoisson, params->maxPoisson, 0.01f, params->poisson);

    if (!params->currentKelvinlet) {
        params->currentKelvinlet = new GrabKelvinlet();
    }
    auto pinchRadio = (new RadioButtonElement("Pinch"))->setChecked(typeid(*(params->currentKelvinlet)) == typeid(PinchKelvinlet));
    auto twistRadio = (new RadioButtonElement("Twist"))->setChecked(typeid(*(params->currentKelvinlet)) == typeid(TwistKelvinlet));
    auto grabRadio = (new RadioButtonElement("Grab"))->setChecked(typeid(*(params->currentKelvinlet)) == typeid(GrabKelvinlet));
    auto scaleRadio = (new RadioButtonElement("Scale"))->setChecked(typeid(*(params->currentKelvinlet)) == typeid(ScaleKelvinlet));

    auto displayAsArrowsCheckbox = (new RadioButtonElement("Arrows"))->setChecked(dataModel->vectorData.displayParameters.displayMode == DisplayedVectorFieldParameters::ARROWS);
    auto displayAsFlowLinesCheckbox = (new RadioButtonElement("Flow lines"))->setChecked(dataModel->vectorData.displayParameters.displayMode == DisplayedVectorFieldParameters::FLOWLINES);
    auto displayNoDisplayCheckbox = (new RadioButtonElement("No display"))->setChecked(dataModel->vectorData.displayParameters.displayMode == DisplayedVectorFieldParameters::NONE);


    auto kelvinletsHistory = new HierarchicalListUI();
    updateKelvinletList(kelvinletsHistory, params);

    // auto editKelvinletCheck = (new CheckboxElement("Edit this Kelvinlet?"))->setChecked(false);
    auto deleteKelvinletButton = new ButtonElement("Delete");

    UI->add(std::vector<UIElement*>{
        radialScaleSlider,
        muSlider,
        poissonSlider,
        createHorizontalGroupUI(std::vector<UIElement*>{pinchRadio, twistRadio, grabRadio, scaleRadio}),
        kelvinletsHistory,
        deleteKelvinletButton,
        createHorizontalGroupUI(std::vector<UIElement*>{displayAsArrowsCheckbox, displayAsFlowLinesCheckbox, displayNoDisplayCheckbox})
    });

    radialScaleSlider->setOnValueChanged([=](float) {
        PainterToolsUI::updateCurrentChartViewWithCurrentKelvinlets(plotter, chartView, dataModel, params, Vector3::invalid, true);
    });
    muSlider->setOnValueChanged([=](float) {
        PainterToolsUI::updateCurrentChartViewWithCurrentKelvinlets(plotter, chartView, dataModel, params, Vector3::invalid, true);
    });
    poissonSlider->setOnValueChanged([=](float) {
        PainterToolsUI::updateCurrentChartViewWithCurrentKelvinlets(plotter, chartView, dataModel, params, Vector3::invalid, true);
    });

    pinchRadio->setOnChecked([=](bool checked) {
        if (!checked) return;
        if (params->currentKelvinlet)
            delete params->currentKelvinlet;
        params->currentKelvinlet = new PinchKelvinlet();
        PainterToolsUI::updateCurrentKelvinlet(params, Vector3::invalid);
    });
    twistRadio->setOnChecked([=](bool checked) {
        if (!checked) return;
        if (params->currentKelvinlet)
            delete params->currentKelvinlet;
        params->currentKelvinlet = new TwistKelvinlet();
        PainterToolsUI::updateCurrentKelvinlet(params, Vector3::invalid);
    });
    grabRadio->setOnChecked([=](bool checked) {
        if (!checked) return;
        if (params->currentKelvinlet)
            delete params->currentKelvinlet;
        params->currentKelvinlet = new GrabKelvinlet();
        PainterToolsUI::updateCurrentKelvinlet(params, Vector3::invalid);
    });
    scaleRadio->setOnChecked([=](bool checked) {
        if (!checked) return;
        if (params->currentKelvinlet)
            delete params->currentKelvinlet;
        params->currentKelvinlet = new ScaleKelvinlet();
        PainterToolsUI::updateCurrentKelvinlet(params, Vector3::invalid);
    });

    displayAsArrowsCheckbox->setOnChecked([=](bool checked) {
        if (!checked) return;
        dataModel->vectorData.displayParameters.displayMode = DisplayedVectorFieldParameters::ARROWS;
        PainterToolsUI::updateCurrentChartViewWithCurrentKelvinlets(plotter, chartView, dataModel, params, Vector3::invalid);
    });
    displayAsFlowLinesCheckbox->setOnChecked([=](bool checked) {
        if (!checked) return;
        dataModel->vectorData.displayParameters.displayMode = DisplayedVectorFieldParameters::FLOWLINES;
        PainterToolsUI::updateCurrentChartViewWithCurrentKelvinlets(plotter, chartView, dataModel, params, Vector3::invalid);
    });
    displayNoDisplayCheckbox->setOnChecked([=](bool checked) {
        if (!checked) return;
        dataModel->vectorData.displayParameters.displayMode = DisplayedVectorFieldParameters::NONE;
        PainterToolsUI::updateCurrentChartViewWithCurrentKelvinlets(plotter, chartView, dataModel, params, Vector3::invalid);
    });


    kelvinletsHistory->setOnItemSelectionChanged([=]() {
        if (kelvinletsHistory->selectedItems().empty()) {
            if (pinchRadio->checked()) params->currentKelvinlet = new PinchKelvinlet();
            else if (twistRadio->checked()) params->currentKelvinlet = new TwistKelvinlet();
            else if (scaleRadio->checked()) params->currentKelvinlet = new ScaleKelvinlet();
            else if (grabRadio->checked()) params->currentKelvinlet = new GrabKelvinlet();
            else std::cerr << "Error: no state for new Kelvinlet" << std::endl;
        }
        else {
            auto item = dynamic_cast<HierarchicalListWidgetItem<Kelvinlet*>*>(kelvinletsHistory->selectedItems()[0]);
            const auto& k = item->getItem();
            params->currentKelvinlet = k;

            params->radialScale = k->radialScale;
            params->mu = k->mu;
            params->poisson = k->v;

            radialScaleSlider->block(); radialScaleSlider->setValue(k->radialScale); radialScaleSlider->unblock();
            muSlider->block(); muSlider->setValue(k->mu); muSlider->unblock();
            poissonSlider->block(); poissonSlider->setValue(k->v); poissonSlider->unblock();

            if (auto asKelvinletPoint = dynamic_cast<KelvinletPoint*>(params->currentKelvinlet)) {
                if (auto asPinch = dynamic_cast<PinchKelvinlet*>(asKelvinletPoint)) {
                    pinchRadio->block(); pinchRadio->setChecked(true); pinchRadio->unblock();
                } else if (auto asTwist = dynamic_cast<TwistKelvinlet*>(asKelvinletPoint)) {
                    twistRadio->block(); twistRadio->setChecked(true); twistRadio->unblock();
                } else if (auto asGrab = dynamic_cast<GrabKelvinlet*>(asKelvinletPoint)) {
                    grabRadio->block(); grabRadio->setChecked(true); grabRadio->unblock();
                } else if (auto asScale = dynamic_cast<ScaleKelvinlet*>(asKelvinletPoint)) {
                    scaleRadio->block(); scaleRadio->setChecked(true); scaleRadio->unblock();
                }

            }
            updateKelvinletList(kelvinletsHistory, params);
        }
    });

    deleteKelvinletButton->setOnClick([=]() {
        auto items = kelvinletsHistory->selectedItems();
        for (auto& _selection : items) {
            auto selection = dynamic_cast<HierarchicalListWidgetItem<Kelvinlet*>*>(_selection);
            auto found = std::find(params->kelvinlets.begin(), params->kelvinlets.end(), selection->itemValue);
            if (found != params->kelvinlets.end()) {
                params->kelvinlets.erase(found);
            }
            updateKelvinletList(kelvinletsHistory, params);            
        }
        PainterToolsUI::updateCurrentChartViewWithCurrentKelvinlets(plotter, chartView, dataModel, params, Vector3::invalid);
    });

    QObject::connect(chartView, &ChartView::clickedOnValue, chartView, [=](const Vector3& pos, bool leftClick, bool rightClick) {
        if (auto k = dynamic_cast<KelvinletPoint*>(params->currentKelvinlet)) {
            if (pos.isValid()) {
                if (rightClick) {
                    params->kelvinletPosition = Vector3::invalid;
                } else {
                    params->kelvinletPosition = pos * dataModel->vectorData.getField().getDimensions();
                }
                k->pos = params->kelvinletPosition;
                k->mu = params->mu;
                k->v = params->poisson;
                k->radialScale = params->radialScale;
            }

            if (!pos.isValid() && k->valid() && !isIn(static_cast<Kelvinlet*>(k), params->kelvinlets)) {
                auto clone = k->clone();
                params->kelvinlets.push_back(clone);
                kelvinletsHistory->addItem(new HierarchicalListWidgetItem<Kelvinlet*>(clone->getShortName(), clone));
                k->reset();
            }
        }
    });

    QObject::connect(chartView, &ChartView::mouseMoved, chartView, [=](const Vector3& pos, const Vector3& prevPos, QMouseEvent* e) {
        PainterToolsUI::updateCurrentChartViewWithCurrentKelvinlets(plotter, chartView, dataModel, params, pos, e->buttons().testFlag(Qt::LeftButton));
    });

    return UI;
}

GridV3& PainterToolsUI::paintKelvinlet(GridV3 &src, const Vector3 &pos, KelvinletToolParams *params)
{
    if (params->currentKelvinlet && params->currentKelvinlet->valid()) {
        src.iterateParallel([&](const Vector3& p) {
            src[p] += params->currentKelvinlet->evaluate(p);
        });
    }
    return src;
}


std::pair<GridV3, GridF> PainterToolsUI::getKelvinletParametersImage(GridV3& img, GridF& alpha, const Vector3& imgScale, Kelvinlet* kelvinlet, const Vector3& pos) {
    if (auto k = dynamic_cast<KelvinletPoint*>(kelvinlet)) {
        Vector3 kCenter = pos.xy();
        img.iterateParallel([&](const Vector3& _p) {
            Vector3 p = _p / imgScale;
            float d = (p - kCenter).norm2();
            if (std::abs(d - (k->radialScale * k->radialScale)) <= (4.f)) {
                img[_p] = Vector3(1, 0, 0);
                alpha[_p] = 1.f;
            } else if (d - (4.f) <= 0) {
                img[_p] = Vector3(1, 0, 0);
                alpha[_p] = 1.f;
            }
        });

        if (k->valid()) {
            Vector3 force;
            if (auto asPinch = dynamic_cast<PinchKelvinlet*>(k)) force = asPinch->force;
            else if (auto asGrab = dynamic_cast<GrabKelvinlet*>(k)) force = asGrab->force;
            else if (auto asTwist= dynamic_cast<TwistKelvinlet*>(k)) force = Vector3(asTwist->force.z(), asTwist->force.z(), 0);
            else if (auto asScale = dynamic_cast<ScaleKelvinlet*>(k)) force = Vector3(asScale->force, asScale->force, 0);

            img = PlotVectorData::drawLine(img, Vector3(0, 0, 1), kCenter.xy() * imgScale, (kCenter + force).xy() * imgScale);
            alpha = PlotVectorData::drawLine(alpha, 1.f, kCenter.xy() * imgScale, (kCenter + force).xy() * imgScale);
        }
    }
    return {img, alpha};
}

void PainterToolsUI::updateCurrentChartViewWithCurrentKelvinlets(AbstractPlotter* plotter, ChartView* chartView, PlotModel* dataModel, KelvinletToolParams* params, const Vector3& mouseRelPos, bool updateCurrentKelvinlet)
{
    if (updateCurrentKelvinlet)
        PainterToolsUI::updateCurrentKelvinlet(params, Vector3::invalid);
    if (!params->displayResultingField) return;
    Vector3i imgSize = Vector3i(300, 300, 1);
    GridV3 img(imgSize);
    GridF alpha(img.getDimensions(), 1.f);
    std::tie(img, alpha) = PlotVectorData::createFieldImageAndAlpha(params->getVectorField(true), img.getDimensions(), Vector3i(20, 20, 1), dataModel->vectorData.displayParameters);

    if (mouseRelPos.isValid()) {
        Vector3 pos = mouseRelPos * dataModel->vectorData.getField().getDimensions();
        if (updateCurrentKelvinlet)
            PainterToolsUI::updateCurrentKelvinlet(params, pos);
    }

    if (auto k = dynamic_cast<KelvinletPoint*>(params->currentKelvinlet))
        PainterToolsUI::getKelvinletParametersImage(img, alpha, (Vector3)imgSize / (Vector3)dataModel->vectorData.field.getDimensions(), k, (k->pos.isValid() ? k->pos : mouseRelPos * dataModel->vectorData.getField().getDimensions()));

    chartView->setOverlay(img, "Kelvinlet placement", alpha, 1000);
    // const GridV3 initial = dataModel->vectorData.field;
    dataModel->vectorData.field = GridV3();
    chartView->setPlotModel(dataModel);
    chartView->update();
    dataModel->vectorData.field = params->getVectorField(true);
    Q_EMIT plotter->updated();
}

Kelvinlet *PainterToolsUI::updateCurrentKelvinlet(KelvinletToolParams *params, const Vector3 &mousePos)
{
    auto& k = params->currentKelvinlet;

    if (mousePos.isValid()) {
        if (auto asPoint = dynamic_cast<KelvinletPoint*>(k)) {
            Vector3 force = mousePos - asPoint->pos;
            if (auto asPinch = dynamic_cast<PinchKelvinlet*>(asPoint)) {
                asPinch->force = force;
            } else if (auto asTwist = dynamic_cast<TwistKelvinlet*>(asPoint)) {
                asTwist->force = Vector3(0, 0, sign(force.x()) * force.norm());
            } else if (auto asGrab = dynamic_cast<GrabKelvinlet*>(asPoint)) {
                asGrab->force = force;
            } else if (auto asScale = dynamic_cast<ScaleKelvinlet*>(asPoint)) {
                asScale->force = sign(force.x()) * force.norm();
            }
        }
    }

    k->mu = params->mu;
    k->radialScale = params->radialScale;
    k->v = params->poisson;
    return k;
}
/*
Kelvinlet *PainterToolsUI::updateSelectedKelvinlet(KelvinletToolParams *params, Kelvinlet* k, const Vector3 &mousePos)
{
    if (mousePos.isValid()) {
        if (auto asPoint = dynamic_cast<KelvinletPoint*>(k)) {
            Vector3 force = mousePos - asPoint->pos;
            if (auto asPinch = dynamic_cast<PinchKelvinlet*>(asPoint)) {
                asPinch->force = force;
            } else if (auto asTwist = dynamic_cast<TwistKelvinlet*>(asPoint)) {
                asTwist->force = Vector3(0, 0, sign(force.x()) * force.norm());
            } else if (auto asGrab = dynamic_cast<GrabKelvinlet*>(asPoint)) {
                asGrab->force = force;
            } else if (auto asScale = dynamic_cast<ScaleKelvinlet*>(asPoint)) {
                asScale->force = sign(force.x()) * force.norm();
            }
        }
    }

    k->mu = params->mu;
    k->radialScale = params->radialScale;
    k->v = params->poisson;
    return k;
}
*/

GridV3 KelvinletToolParams::getVectorField(bool takeIntoAccountCurrentKelvinlet) const
{
    GridV3 resultingVectorField = this->getInitialVectorField();
    auto evaluatingKelvinlets = std::vector<Kelvinlet*>();
    for (const auto& k : this->kelvinlets) {
        if (k && k->valid())
            evaluatingKelvinlets.push_back(k);
    }
    for (const auto& k : this->additional_kelvinlets) {
        if (k && k->valid())
            evaluatingKelvinlets.push_back(k);
    }
    if (takeIntoAccountCurrentKelvinlet && this->currentKelvinlet->valid() && !isIn(static_cast<Kelvinlet*>(this->currentKelvinlet), evaluatingKelvinlets))
        evaluatingKelvinlets.push_back(this->currentKelvinlet->clone());

    resultingVectorField.iterateParallel([&](const Vector3i& p) {
        for (const auto& k : evaluatingKelvinlets)
            resultingVectorField[p] += k->evaluate(p);
    });
    return resultingVectorField;
}
