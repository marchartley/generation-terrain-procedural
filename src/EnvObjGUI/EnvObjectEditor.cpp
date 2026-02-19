#include "EnvObjectEditor.h"

#include "GUIElements/ImageViewerOptionsUI.h"

#include "serialization/Serializer.h"

EnvObjectEditor::EnvObjectEditor(const std::string& name, QWidget* parent) : EnvObjectEditor(name, new ChartView(new Chart()), parent)
{}
EnvObjectEditor::EnvObjectEditor(const std::string& name, ChartView* chartView, QWidget* parent)
    : ImageViewer(name, chartView, parent)
{
    dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::red, Vector3::white, Vector3::green});

    kelvinletParams.displayResultingField = false;

    bodyParameters.pinchK = new PinchKelvinletCurve();
    bodyParameters.twistK = new TwistKelvinletCurve();
    bodyParameters.grabK = new GrabKelvinletCurve();
    bodyParameters.scaleK = new ScaleKelvinletCurve();

    kelvinletParams.additional_kelvinlets.push_back(bodyParameters.pinchK);
    kelvinletParams.additional_kelvinlets.push_back(bodyParameters.twistK);
    kelvinletParams.additional_kelvinlets.push_back(bodyParameters.grabK);
    kelvinletParams.additional_kelvinlets.push_back(bodyParameters.scaleK);
}

EnvObjectEditor *EnvObjectEditor::updateToolsInterface()
{
    this->toolsInterface->clear();
    this->toolsInterface->add(PainterToolsUI::createKelvinletToolsUI(this, this->chartView, this->dataModel, &this->kelvinletParams));

    InterfaceUI* bodyKelvinletUI = new InterfaceUI(new QVBoxLayout(), false, "Kelvinlet body editor");

    auto pinchForceSlider = new SliderElement("Pinch force", bodyParameters.minForce, bodyParameters.maxForce, .01f, bodyParameters.pinchK->force);
    auto twistForceSlider = new SliderElement("Twist force", bodyParameters.minForce, bodyParameters.maxForce, .01f, bodyParameters.twistK->force);
    auto grabForceSlider = new SliderElement("Grab force", bodyParameters.minForce, bodyParameters.maxForce, .01f, bodyParameters.grabK->force);
    auto scaleForceSlider = new SliderElement("Scale force", bodyParameters.minForce, bodyParameters.maxForce, .01f, bodyParameters.scaleK->force);

    auto simulationModeCheckbox = new CheckboxElement("Simulation", depositionSimulationDisplay);

    auto angleInitialFlow = new AngleElement("Flow angle");
    auto strengthInitialFlow = (new SliderElement("Flow strength", 0.f, 3.f, .01f))->setValue(1.f);

    auto validationButton = new ButtonElement("Validate");

    if (this->currentObject->isCurve() || this->currentObject->isArea()) {
        bodyKelvinletUI->add(std::vector<UIElement*>{
            pinchForceSlider,
            twistForceSlider,
            grabForceSlider,
            scaleForceSlider
        });
    }
    bodyKelvinletUI->add(std::vector<UIElement*>{
        simulationModeCheckbox,
        angleInitialFlow,
        strengthInitialFlow,
        validationButton
    });

    this->toolsInterface->add(bodyKelvinletUI);




    pinchForceSlider->setOnValueChanged([=](float newValue) {
        this->updateCurrentChartViewWithCurrentKelvinlets(Vector3::invalid, false);
    });

    twistForceSlider->setOnValueChanged([=](float newValue) {
        this->updateCurrentChartViewWithCurrentKelvinlets(Vector3::invalid, false);
    });

    grabForceSlider->setOnValueChanged([=](float newValue) {
        this->updateCurrentChartViewWithCurrentKelvinlets(Vector3::invalid, false);
    });

    scaleForceSlider->setOnValueChanged([=](float newValue) {
        this->updateCurrentChartViewWithCurrentKelvinlets(Vector3::invalid, false);
    });

    QObject::connect(this, &AbstractPlotter::updated, this, [=]() {
        if (this->depositionSimulationDisplay) {
            depositionGrid = GridF(100, 100, 1);
            this->simulateDeposition(depositionGrid);
            this->setOverlay(depositionGrid, GridF(), "deposition simulation", -1);
            this->showOverlay("deposition simulation");
        } else {
            this->hideOverlay("deposition simulation");
        }
    });

    angleInitialFlow->setOnValueChanged([=](float newAngle) {
        this->kelvinletParams.temporaryVectorData.field = GridV3(this->kelvinletParams.temporaryVectorData.field.getDimensions(), Vector3(0, strengthInitialFlow->value(), 0).rotated(Vector3(0, 0, -deg2rad(newAngle))));
        this->updateCurrentChartViewWithCurrentKelvinlets(Vector3::invalid, false);
    });
    strengthInitialFlow->setOnValueChanged([=](float newStrength) {
        this->kelvinletParams.temporaryVectorData.field = GridV3(this->kelvinletParams.temporaryVectorData.field.getDimensions(), Vector3(0, newStrength, 0).rotated(Vector3(0, 0, -deg2rad(angleInitialFlow->value()))));
        this->updateCurrentChartViewWithCurrentKelvinlets(Vector3::invalid, false);
    });

    validationButton->setOnPressed([=]() {
        validateEnvObject();
    });


    QObject::connect(chartView, &ChartView::mouseMoved, chartView, [=](const Vector3& pos, const Vector3& prevPos, QMouseEvent* e) {
        this->updateCurrentChartViewWithCurrentKelvinlets(pos, e->buttons().testFlag(Qt::LeftButton));
    });



    return this;
}

EnvObjectEditor *EnvObjectEditor::addEnvObject(EnvObject *envObj)
{
    this->currentObject = envObj->clone();

    if (auto asPoint = dynamic_cast<EnvPoint*>(this->currentObject)) {
        asPoint->position = Vector3(50, 50, 0);
    } else if (auto asCurve = dynamic_cast<EnvCurve*>(this->currentObject)) {
        asCurve->curve = BSpline({Vector3(10, 10), Vector3(60, 30), Vector3(30, 60), Vector3(90, 90)});

        bodyParameters.pinchK->curve = asCurve->curve;
        bodyParameters.twistK->curve = asCurve->curve;
        bodyParameters.grabK->curve = asCurve->curve;
        bodyParameters.scaleK->curve = asCurve->curve;

        bodyParameters.pinchK->radialScale = asCurve->width;
        bodyParameters.twistK->radialScale = asCurve->width;
        bodyParameters.grabK->radialScale = asCurve->width;
        bodyParameters.scaleK->radialScale = asCurve->width;
    } else if (auto asArea = dynamic_cast<EnvArea*>(this->currentObject)) {
        asArea->curve = ShapeCurve({Vector3(30, 30), Vector3(30, 70), Vector3(70, 70), Vector3(60, 40)});

        bodyParameters.pinchK->curve = asArea->curve;
        bodyParameters.twistK->curve = asArea->curve;
        bodyParameters.grabK->curve = asArea->curve;
        bodyParameters.scaleK->curve = asArea->curve;

        bodyParameters.pinchK->radialScale = asArea->width;
        bodyParameters.twistK->radialScale = asArea->width;
        bodyParameters.grabK->radialScale = asArea->width;
        bodyParameters.scaleK->radialScale = asArea->width;
    }


    this->addImage(displayEnvObject().first);
    return this;
}

std::pair<GridV3, GridF> EnvObjectEditor::displayEnvObject() const
{
    Vector3 imgScale = Vector3(2.f, 2.f, 1.f);
    GridV3 img(100 * imgScale.x(), 100 * imgScale.y(), 1, this->dataModel->vectorData.displayParameters.backgroundColor);
    GridF alpha(img.getDimensions());
    if (auto asPoint = dynamic_cast<EnvPoint*>(this->currentObject)) {
        float crossLength = 3;
        Vector3 pos = asPoint->position * imgScale;

        PlotVectorData::drawLine(img, Vector3(0, 0, 1), pos - Vector3(crossLength / 2, crossLength / 2) * imgScale, pos + Vector3(crossLength / 2, crossLength / 2) * imgScale);
        PlotVectorData::drawLine(img, Vector3(0, 0, 1), pos - Vector3(crossLength / 2, -crossLength / 2) * imgScale, pos + Vector3(crossLength / 2, -crossLength / 2) * imgScale);

        PlotVectorData::drawCircle(img, Vector3(1, 0, 0), pos, asPoint->radius);

        PlotVectorData::drawLine(alpha, 1.f, pos - Vector3(crossLength / 2, crossLength / 2) * imgScale, pos + Vector3(crossLength / 2, crossLength / 2) * imgScale);
        PlotVectorData::drawLine(alpha, 1.f, pos - Vector3(crossLength / 2, -crossLength / 2) * imgScale, pos + Vector3(crossLength / 2, -crossLength / 2) * imgScale);

        PlotVectorData::drawCircle(alpha, 1.f, pos, asPoint->radius);

    } else if (auto asCurve = dynamic_cast<EnvCurve*>(this->currentObject)) {
        BSpline curve = asCurve->curve;
        img.iterateParallel([&](const Vector3& p) {
            float distToCurve = curve.estimateDistanceFrom(p / imgScale);

            if (distToCurve < 1.f) {
                img(p) = Vector3(0, 0, 1);
                alpha(p) = 1.f;
            }
            if (abs(distToCurve - asCurve->width) < 1.f) {
                img(p) = Vector3(1, 0, 0);
                alpha(p) = 1.f;
            }
        });
    } else if (auto asArea = dynamic_cast<EnvArea*>(this->currentObject)) {
        ShapeCurve curve = asArea->curve;
        img.iterateParallel([&](const Vector3& p) {
            float distToCurve = curve.estimateDistanceFrom(p / imgScale);

            if (abs(distToCurve) < 1.f) {
                img(p) = Vector3(0, 0, 1);
                alpha(p) = 1.f;
            }
        });
    }
    return {img, alpha};
}

GridF& EnvObjectEditor::simulateDeposition(GridF &currentState, int iterations)
{
    const GridV3& flow = this->dataModel->vectorData.getField();
    currentState = currentState.resize(flow.getDimensions());
    currentState.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::WRAPPED_VALUE;

    for (int i = 0; i < iterations; i++) {
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


GridV3 EnvObjectEditor::getVectorFieldWithRotation(bool takeIntoAccountCurrentKelvinlet) const
{
    this->validateEnvObject();
    GridV3 resultingVectorField = this->kelvinletParams.getInitialVectorField();
    currentObject->computeFlowModification(resultingVectorField);
    return resultingVectorField;
}

GridV3 EnvObjectEditor::getVectorFieldWithRotationForEnvPoint(bool takeIntoAccountCurrentKelvinlet) const
{
    GridV3 resultingVectorField = this->kelvinletParams.getInitialVectorField();
    auto evaluatingRelativeKelvinlets = std::vector<RelativeKelvinlet>();
    auto evaluatingKelvinlets = std::set<Kelvinlet*>();
    for (const auto& k : this->kelvinletParams.kelvinlets) {
        if (k && k->valid()) {
            evaluatingRelativeKelvinlets.push_back(RelativeKelvinlet(k, dynamic_cast<EnvPoint*>(currentObject)->position));
            evaluatingKelvinlets.insert(k);
        }
    }
    for (const auto& k : this->kelvinletParams.additional_kelvinlets) {
        if (k && k->valid()) {
            evaluatingRelativeKelvinlets.push_back(RelativeKelvinlet(k, dynamic_cast<EnvPoint*>(currentObject)->position));
            evaluatingKelvinlets.insert(k);
        }
    }
    if (takeIntoAccountCurrentKelvinlet && this->kelvinletParams.currentKelvinlet->valid() && !isIn(this->kelvinletParams.currentKelvinlet, evaluatingKelvinlets))
        evaluatingRelativeKelvinlets.push_back(RelativeKelvinlet(this->kelvinletParams.currentKelvinlet->clone(), dynamic_cast<EnvPoint*>(currentObject)->position));

    float flowAngle = resultingVectorField.at(dynamic_cast<EnvPoint*>(currentObject)->position).getSignedAngleWith(Vector3(1, 0, 0));
    float flowStrength = resultingVectorField.at(dynamic_cast<EnvPoint*>(currentObject)->position).norm();
    resultingVectorField.iterateParallel([&](const Vector3i& p) {
        for (const auto& k : evaluatingRelativeKelvinlets)
            resultingVectorField[p] += k.evaluate(p, flowAngle, flowStrength);
    });
    return resultingVectorField;
}

GridV3 EnvObjectEditor::getVectorFieldWithRotationForEnvCurve(bool takeIntoAccountCurrentKelvinlet) const
{
    GridV3 resultingVectorField = this->kelvinletParams.getInitialVectorField();
    return resultingVectorField;
}

GridV3 EnvObjectEditor::getVectorFieldWithRotationForEnvArea(bool takeIntoAccountCurrentKelvinlet) const
{
    GridV3 resultingVectorField = this->kelvinletParams.getInitialVectorField();
    return resultingVectorField;
}

void EnvObjectEditor::updateCurrentChartViewWithCurrentKelvinlets(const Vector3& mouseRelPos, bool updateCurrentKelvinlet)
{
    if (updateCurrentKelvinlet)
        PainterToolsUI::updateCurrentKelvinlet(&this->kelvinletParams, Vector3::invalid);
    Vector3i imgSize = Vector3i(300, 300, 1);
    GridV3 img(imgSize);
    GridF alpha(img.getDimensions(), 1.f);
    std::tie(img, alpha) = PlotVectorData::createFieldImageAndAlpha(this->getVectorFieldWithRotation(true), img.getDimensions(), Vector3i(20, 20, 1), dataModel->vectorData.displayParameters);

    if (mouseRelPos.isValid()) {
        Vector3 pos = mouseRelPos * dataModel->vectorData.getField().getDimensions();
        if (updateCurrentKelvinlet)
            PainterToolsUI::updateCurrentKelvinlet(&this->kelvinletParams, pos);
    }

    if (auto k = dynamic_cast<KelvinletPoint*>(this->kelvinletParams.currentKelvinlet))
        PainterToolsUI::getKelvinletParametersImage(img, alpha, (Vector3)imgSize / (Vector3)dataModel->vectorData.field.getDimensions(), k, (k->pos.isValid() ? k->pos : mouseRelPos * dataModel->vectorData.getField().getDimensions()));

    chartView->setOverlay(img, "Kelvinlet placement", alpha, 1000);
    // const GridV3 initial = dataModel->vectorData.field;
    dataModel->vectorData.field = GridV3();
    chartView->setPlotModel(dataModel);
    chartView->update();
    dataModel->vectorData.field = this->getVectorFieldWithRotation(true);
    Q_EMIT this->updated();
}

EnvObject *EnvObjectEditor::validateEnvObject() const
{
    if (auto asPoint = dynamic_cast<EnvPoint*>(currentObject)) {
        // asPoint->mainKelvinlets = this->kelvinletParams.kelvinlets;
        asPoint->mainKelvinlets.clear();
        for (auto& k : kelvinletParams.kelvinlets) {
            asPoint->mainKelvinlets.push_back(k->clone()->translate(-asPoint->position));
        }
    }
    else if (auto asCurve = dynamic_cast<EnvCurve*>(currentObject)) {
        // asCurve->startingPointKelvinlets = ...;
        // asCurve->endingPointKelvinlets = ...;
        asCurve->curveKelvinlets = this->kelvinletParams.additional_kelvinlets;
    }
    else if (auto asArea = dynamic_cast<EnvArea*>(currentObject)) {
        asArea->curveKelvinlets = this->kelvinletParams.additional_kelvinlets;
    }

    // nlohmann::json json = currentObject;
    // std::cout << "My object in JSON form: \n" << json.dump(1) << std::endl;

    return currentObject;
}
