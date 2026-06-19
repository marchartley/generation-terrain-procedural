#include "EnvObjectEditor.h"

#include "GUIElements/ImageViewerOptionsUI.h"

#include "serialization/Serializer.h"

EnvObjectEditor::EnvObjectEditor(const std::string& name, QWidget* parent)
    : ImageViewer(name, parent)
{
    dataModel->imageData.displayParameters.colorRamp = CatmullRomSpline({Vector3::red, Vector3::white, Vector3::green});

    kelvinletParams.displayResultingField = false;

    /*QTimer* animationTimer = new QTimer();
    animationTimer->connect(animationTimer, &QTimer::timeout, this, [=]() {
    });
    animationTimer->setInterval(20);
    animationTimer->start();*/

    this->animate([=]() {
        if (this->currentObject) {
            displayDepositionSimulation();
            if (this->animating) {
                updateCurrentChartViewWithCurrentKelvinlets(Vector3::invalid, false);
            }
        }
        return true;
    });
}

EnvObjectEditor& EnvObjectEditor::updateToolsInterface()
{
    auto updateView = [=](const Vector3& mousePos, bool updateCurrentKelvinlet) {
        this->updateCurrentChartViewWithCurrentKelvinlets(mousePos, updateCurrentKelvinlet);
    };
    auto updateField = [=](bool useCurrentKelvinlet) {
        validateEnvObject(useCurrentKelvinlet);
        return this->getVectorFieldWithRotation(useCurrentKelvinlet);
    };
    this->toolsInterface->clear();
    this->toolsInterface->add(std::move(PainterToolsUI::createKelvinletToolsUI(this, &this->kelvinletParams, updateView, updateField)));

    auto bodyKelvinletUI = new InterfaceUI(InterfaceUI::VERTICAL, false, "Kelvinlet body editor");

    auto simulationModeCheckbox = new CheckboxElement("Simulation", depositionSimulationDisplay);
    auto animatedModeCheckbox = new CheckboxElement("Animate", animating);

    auto anchorSelectionCombobox = new ComboboxElement("Anchor");

    auto angleInitialFlow = new AngleElement("Flow angle");
    auto strengthInitialFlow = new SliderElement("Flow strength", 0.f, 3.f, .01f);

    if (this->currentObject) {
        if (this->currentObject->getDefinition()->isPoint()) {
            anchorSelectionCombobox->addChoice(new ComboboxLineElement<KELVINLET_ANCHOR_POINT>("Center", MAIN), true);
        }
        else if (this->currentObject->getDefinition()->isCurve()) {
            anchorSelectionCombobox->addChoice(new ComboboxLineElement<KELVINLET_ANCHOR_POINT>("Start", START), true);
            anchorSelectionCombobox->addChoice(new ComboboxLineElement<KELVINLET_ANCHOR_POINT>("End", END));
            anchorSelectionCombobox->addChoice(new ComboboxLineElement<KELVINLET_ANCHOR_POINT>("Curve", CURVE));
        } else {
            anchorSelectionCombobox->addChoice(new ComboboxLineElement<KELVINLET_ANCHOR_POINT>("Curve", CURVE), true);
        }
        this->currentAnchorPoint = dynamic_cast<ComboboxLineElement<KELVINLET_ANCHOR_POINT>*>(anchorSelectionCombobox->choices[0])->value;
        this->kelvinletAnchors[kelvinletParams.currentKelvinlet] = {currentAnchorPoint, Vector3::invalid};
    }

    angleInitialFlow->setAngle(0.f);
    strengthInitialFlow->setValue(1.f);

    auto validationButton = new ButtonElement("Validate", [=]() { validateEnvObject(false, true, true); });

    bodyKelvinletUI->add({
        anchorSelectionCombobox,
        createHorizontalGroup(std::vector<UIElement*>{simulationModeCheckbox, animatedModeCheckbox}),
        angleInitialFlow,
        strengthInitialFlow,
        validationButton,
    });

    this->toolsInterface->add(std::move(bodyKelvinletUI));

    angleInitialFlow->setOnValueChanged([=](float newAngle) {
        this->kelvinletParams.temporaryVectorData.field = GridV3(this->kelvinletParams.temporaryVectorData.field.getDimensions(), Vector3(0, strengthInitialFlow->value(), 0).rotated(Vector3(0, 0, -deg2rad(newAngle))));
        this->updateCurrentChartViewWithCurrentKelvinlets(Vector3::invalid, false);
    });
    strengthInitialFlow->setOnValueChanged([=](float newStrength) {
        this->kelvinletParams.temporaryVectorData.field = GridV3(this->kelvinletParams.temporaryVectorData.field.getDimensions(), Vector3(0, newStrength, 0).rotated(Vector3(0, 0, -deg2rad(angleInitialFlow->value()))));
        this->updateCurrentChartViewWithCurrentKelvinlets(Vector3::invalid, false);
    });

    kelvinletParams.setOnNewKelvinlet([=](Kelvinlet* k) {
        auto anchorType = currentAnchorPoint;
        switch (anchorType) {
            case KELVINLET_ANCHOR_POINT::MAIN :
                this->kelvinletAnchors[k] = std::make_pair(anchorType, dynamic_cast<EnvPointInstance*>(currentObject)->position);
                break;
            case KELVINLET_ANCHOR_POINT::START :
                this->kelvinletAnchors[k] = std::make_pair(anchorType, dynamic_cast<EnvCurveInstance*>(currentObject)->curve.front());
                break;
            case KELVINLET_ANCHOR_POINT::END :
                this->kelvinletAnchors[k] = std::make_pair(anchorType, dynamic_cast<EnvCurveInstance*>(currentObject)->curve.back());
                break;
            case KELVINLET_ANCHOR_POINT::CURVE :
                this->kelvinletAnchors[k] = std::make_pair(anchorType, Vector3::origin);
                break;
            default:
                break;
        }
        this->validateEnvObject(false, true, false);
    });
    kelvinletParams.setOnDeletedKelvinlet([=](Kelvinlet* k) {
        this->validateEnvObject(false, true, false);
    });
    kelvinletParams.setOnModifiedKelvinlet([=](Kelvinlet* k) {
        this->validateEnvObject(false, true, false);
    });

    anchorSelectionCombobox->setOnSelectionChanged([=] (int) {
        this->currentAnchorPoint = anchorSelectionCombobox->getSelection<KELVINLET_ANCHOR_POINT>();
    });
    return *this;
}

EnvObjectEditor& EnvObjectEditor::addEnvObject(EnvObject *envObj)
{
    this->currentObject = envObj->instantiate();

    bodyParameters.resetKelvinlets();
    kelvinletParams.resetKelvinlets();

    if (auto list = dynamic_cast<HierarchicalListUI*>(this->toolsInterface->findByName("kelvinlets", true))) {
        list->clear();
    }

    if (auto asPoint = dynamic_cast<EnvPointInstance*>(this->currentObject)) {
        asPoint->position = Vector3(50, 50, 0);
        this->objectScale = 100.f / (2.f * asPoint->getDefinition()->radius);

        for (auto& k : asPoint->getDefinition()->mainKelvinlets) {
            auto clone = k->clone();
            clone->scale(objectScale);
            clone->translate(asPoint->position);
            kelvinletParams.kelvinlets.push_back(clone);
            this->kelvinletAnchors[clone] = {MAIN, asPoint->position};
        }
    } else if (auto asCurve = dynamic_cast<EnvCurveInstance*>(this->currentObject)) {
        asCurve->curve = CatmullRomSpline({Vector3(10, 10), Vector3(60, 30), Vector3(30, 60), Vector3(90, 90)});
        this->objectScale = asCurve->curve.length() / asCurve->getDefinition()->length;
    } else if (auto asArea = dynamic_cast<EnvAreaInstance*>(this->currentObject)) {
        asArea->curve = Contour({Vector3(30, 30), Vector3(30, 70), Vector3(70, 70), Vector3(60, 40)});
        this->objectScale = asArea->curve.computeArea() / (asArea->getDefinition()->width * asArea->getDefinition()->length);
    }

    this->addImage(displayEnvObject().first);
    return *this;
}

std::pair<GridV3, GridF> EnvObjectEditor::displayEnvObject() const
{
    Vector3 imgScale = Vector3(2.f, 2.f, 1.f);
    GridV3 img(100 * imgScale.x(), 100 * imgScale.y(), 1, this->dataModel->vectorData.displayParameters.backgroundColor);
    GridF alpha(img.getDimensions());
    if (auto asPoint = dynamic_cast<EnvPointInstance*>(this->currentObject)) {
        float crossLength = asPoint->getDefinition()->radius;
        Vector3 pos = asPoint->position * imgScale;

        PlottingUtils::drawLine(img, Vector3(0, 0, 1), pos - Vector3(crossLength * .5f, crossLength * .5f) * imgScale, pos + Vector3(crossLength * .5f, crossLength * .5f) * imgScale);
        PlottingUtils::drawLine(img, Vector3(0, 0, 1), pos - Vector3(crossLength * .5f, -crossLength * .5f) * imgScale, pos + Vector3(crossLength * .5f, -crossLength * .5f) * imgScale);

        PlottingUtils::drawCircle(img, Vector3(1, 0, 0), pos, asPoint->getDefinition()->radius * this->objectScale);

        PlottingUtils::drawLine(alpha, 1.f, pos - Vector3(crossLength * .5f, crossLength * .5f) * imgScale, pos + Vector3(crossLength * .5f, crossLength * .5f) * imgScale);
        PlottingUtils::drawLine(alpha, 1.f, pos - Vector3(crossLength * .5f, -crossLength * .5f) * imgScale, pos + Vector3(crossLength * .5f, -crossLength * .5f) * imgScale);

        PlottingUtils::drawCircle(alpha, 1.f, pos, asPoint->getDefinition()->radius * this->objectScale);

    } else if (auto asCurve = dynamic_cast<EnvCurveInstance*>(this->currentObject)) {
        CatmullRomSpline curve = asCurve->curve.scaled(imgScale);
        PlottingUtils::drawCurve(img, Vector3(0, 0, 1), curve);
        PlottingUtils::drawCurve(alpha, 1.f, curve);

    } else if (auto asArea = dynamic_cast<EnvAreaInstance*>(this->currentObject)) {
        Contour curve(asArea->curve.curve->clone()->scale(imgScale));
        PlottingUtils::drawCurve(img, Vector3(0, 0, 1), *curve.curve);
        PlottingUtils::drawCurve(alpha, 1.f, *curve.curve);
    }
    return {img, alpha};
}

GridF& EnvObjectEditor::simulateDeposition(GridF &currentState, int iterations)
{
    const GridV3& flow = this->dataModel->vectorData.getField();
    currentState = currentState.resize(flow.getDimensions());
    currentState.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::WRAPPED_VALUE;

    for (int iter = 0; iter < iterations; iter++) {
        currentState.iterateParallel([&](size_t i){ currentState[i] = std::max(currentState[i] + 1.f, 0.f); });

        currentState = currentState.warpWith(flow);
        currentState.iterate([&](size_t x, size_t y, size_t){
            if (x != 0 && x != currentState.sizeX - 1 && y != 0 && y != currentState.sizeY - 1) return;
            currentState(x, y) = std::max({currentState.at(x, y), currentState.at(x-1, y), currentState.at(x+1, y), currentState.at(x, y-1), currentState.at(x, y+1)});
        });
        currentState.iterateParallel([&](size_t i){ currentState[i] = currentState[i] - 1.f; });
    }
    currentState.iterateParallel([&](size_t i){ currentState[i] = std::max(currentState[i] + .5f, 0.f); });
    return currentState;
}


GridV3 EnvObjectEditor::getVectorFieldWithRotation(bool)
{
    GridV3 resultingVectorField = this->kelvinletParams.getInitialVectorField();
    currentObject->computeFlowModification(resultingVectorField, this->objectScale);
    return resultingVectorField;
}

void EnvObjectEditor::updateCurrentChartViewWithCurrentKelvinlets(const Vector3& mouseRelPos, bool updateCurrentKelvinlet)
{
    this->animateEnvObject(animating);
    if (updateCurrentKelvinlet)
        PainterToolsUI::updateCurrentKelvinlet(&this->kelvinletParams, Vector3::invalid);

    Vector3i imgSize = Vector3i(300, 300, 1);
    GridV3 img(imgSize);
    GridF alpha(img.getDimensions(), 1.f);

    this->validateEnvObject(true, false);
    GridV3 field = this->getVectorFieldWithRotation();

    std::tie(img, alpha) = PlotVectorData::createFieldImageAndAlpha(field, img.getDimensions(), Vector3i(20, 20, 1), dataModel->vectorData.displayParameters);

    if (mouseRelPos.isValid()) {
        Vector3 pos = mouseRelPos * dataModel->vectorData.getField().getDimensions();
        if (updateCurrentKelvinlet)
            PainterToolsUI::updateCurrentKelvinlet(&this->kelvinletParams, pos);
    }

    if (auto k = dynamic_cast<KelvinletPoint*>(this->kelvinletParams.currentKelvinlet)) {
        PainterToolsUI::getKelvinletParametersImage(img, alpha, (Vector3)imgSize / (Vector3)dataModel->vectorData.field.getDimensions(), k, (k->pos.isValid() ? k->pos : mouseRelPos * dataModel->vectorData.getField().getDimensions()));
    }

    chartView->setOverlay(img, "Kelvinlet placement", alpha, 1000);
    dataModel->vectorData.field = GridV3();
    chartView->setPlotModel(dataModel);
    chartView->update();
    dataModel->vectorData.field = field;

    emitUpdate();
}

EnvObject *EnvObjectEditor::validateEnvObject(bool takeIntoAccountCurrentKelvinlet, bool sendSignal, bool sendDefinitiveObject)
{
    std::vector<Kelvinlet*> evaluatedKelvinlets = kelvinletParams.kelvinlets;
    if (takeIntoAccountCurrentKelvinlet && this->kelvinletParams.currentKelvinlet->valid() && !isIn(kelvinletParams.currentKelvinlet, kelvinletParams.kelvinlets)) {
        evaluatedKelvinlets.push_back(this->kelvinletParams.currentKelvinlet);
        auto anchorType = currentAnchorPoint;
        switch (anchorType) {
            case KELVINLET_ANCHOR_POINT::MAIN:
                this->kelvinletAnchors[this->kelvinletParams.currentKelvinlet] = std::make_pair(anchorType, dynamic_cast<EnvPointInstance*>(currentObject)->position);
                break;
            case KELVINLET_ANCHOR_POINT::START:
                this->kelvinletAnchors[this->kelvinletParams.currentKelvinlet] = std::make_pair(anchorType, dynamic_cast<EnvCurveInstance*>(currentObject)->curve.front());
                break;
            case KELVINLET_ANCHOR_POINT::END:
                this->kelvinletAnchors[this->kelvinletParams.currentKelvinlet] = std::make_pair(anchorType, dynamic_cast<EnvCurveInstance*>(currentObject)->curve.back());
                break;
            case KELVINLET_ANCHOR_POINT::CURVE:
                this->kelvinletAnchors[this->kelvinletParams.currentKelvinlet] = std::make_pair(anchorType, Vector3::origin);
                break;
            default:
                break;
        }
    }

    currentObject->getDefinition()->clearKelvinlets();
    if (auto asPoint = dynamic_cast<EnvPoint*>(currentObject->getDefinition())) {
        for (auto& k : evaluatedKelvinlets) {
            auto newK = k->clone();
            newK->translate(-kelvinletAnchors[k].second);
            newK->scale(1.f / this->objectScale);
            asPoint->mainKelvinlets.push_back(newK);
        }
    }
    else if (auto asCurve = dynamic_cast<EnvCurve*>(currentObject->getDefinition())) {
        for (auto& k : evaluatedKelvinlets) {
            auto newK = k->clone();
            newK->translate(-kelvinletAnchors[k].second);
            newK->scale(1.f / this->objectScale);
            if (kelvinletAnchors.count(k) > 0 && kelvinletAnchors.at(k).first == START) {
                // newK->translate(-kelvinletAnchors[k].second);
                asCurve->startingPointKelvinlets.push_back(newK);
            } else if (kelvinletAnchors.count(k) > 0 && kelvinletAnchors.at(k).first == END) {
                asCurve->endingPointKelvinlets.push_back(newK);
            } else if (kelvinletAnchors.count(k) > 0 && kelvinletAnchors.at(k).first == CURVE) {
                if (auto pointKelvinlet = dynamic_cast<KelvinletPoint*>(newK)) {
                    asCurve->curveKelvinlets.push_back(pointKelvinlet->cloneToCurveKelvinlet());
                }
            }
        }
    }
    else if (auto asArea = dynamic_cast<EnvArea*>(currentObject->getDefinition())) {
        for (auto& k : evaluatedKelvinlets) {
            auto newK = k->clone();
            newK->translate(-kelvinletAnchors[k].second);
            newK->scale(1.f / this->objectScale);
            if (auto pointKelvinlet = dynamic_cast<KelvinletPoint*>(newK)) {
                asArea->curveKelvinlets.push_back(pointKelvinlet->cloneToCurveKelvinlet());
            }
        }
    }

    if (sendSignal) {
        emitObjectModified(currentObject->getDefinition());
    }
    if (sendDefinitiveObject) {
        emitObjectSaved(currentObject->getDefinition());
    }
    return currentObject->getDefinition();
}

void EnvObjectEditor::animateEnvObject(bool animate)
{
    if (animate) {
        if (animationFrame == 0) {
            AABBox bounds(Vector3::origin, dataModel->vectorData.field.getDimensions().xy());
            if (dynamic_cast<EnvPointInstance*>(this->currentObject)) {
                verticesTargets = {Vector3::random(bounds)};
            }
            else if (auto asCurve = dynamic_cast<EnvCurveInstance*>(this->currentObject)) {
                verticesTargets.resize(asCurve->curve.size());
                for (auto& p : verticesTargets)
                    p = Vector3::random(bounds);
            }
            else if (auto asArea = dynamic_cast<EnvAreaInstance*>(this->currentObject)) {
                verticesTargets.resize(asArea->curve.curve->size());
                for (auto& p : verticesTargets)
                    p = Vector3::random(bounds);
            }
        }
        float remainingFrames = 50 - animationFrame;
        if (auto asPoint = dynamic_cast<EnvPointInstance*>(this->currentObject)) {
            auto& p = asPoint->position;
            p += (verticesTargets[0] - p) * (interpolation::smooth(remainingFrames / 50.f)) / remainingFrames;
        }
        else if (auto asCurve = dynamic_cast<EnvCurveInstance*>(this->currentObject)) {
            auto newCurve = asCurve->curve;
            for (size_t i = 0; i < newCurve.size(); i++) {
                auto& p = newCurve[i];
                p += (verticesTargets[i] - p) * (interpolation::smooth(remainingFrames / 50.f)) / remainingFrames;
            }
            asCurve->updateCurve(newCurve);
        }
        else if (auto asArea = dynamic_cast<EnvAreaInstance*>(this->currentObject)) {
            auto newCurve = asArea->curve;
            for (size_t i = 0; i < newCurve.curve->size(); i++) {
                auto& p = newCurve.curve->get(i);
                p += (verticesTargets[i] - p) * (interpolation::smooth(remainingFrames / 50.f)) / remainingFrames;
            }
            asArea->updateCurve(*newCurve.curve);
        }
        animationFrame = (animationFrame + 1) % 50;
    }
    this->addImage(displayEnvObject().first);
}

void EnvObjectEditor::displayDepositionSimulation()
{
    if (this->depositionSimulationDisplay) {
        depositionGrid = GridF(100, 100, 1);
        this->simulateDeposition(depositionGrid);
        this->setOverlay(depositionGrid, GridF(), "deposition simulation", -1);
        this->showOverlay("deposition simulation");
    } else {
        this->hideOverlay("deposition simulation");
    }
}

void BodyKelvinletParameters::resetKelvinlets()
{
    for (auto& rK : relativeKelvinlets)
        delete rK;
    relativeKelvinlets.clear();
}
