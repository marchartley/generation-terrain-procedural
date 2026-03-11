#include "SnakeSegmentationEditor.h"

#include "GUIElements/PainterToolsUI.h"
#include "GUIElements/ImageViewerOptionsUI.h"

#include "EnvObject/PositionOptimizer.h"

#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"

SnakeSegmentationEditor::SnakeSegmentationEditor(const std::string& name, QWidget* parent)
    : ImageViewer(name, parent)
{
    this->painterParameters.RGBimage = false;

    this->snakeParameters.snake = SnakeSegmentation();
    this->snakeParameters.params = new SnakeSegmentationParameters();
    this->snakeParameters.field = new SnakeImageFieldImplicit();

    this->setSnakeImage(GridF::perlin(Vector3i(100, 100, 1), Vector3(5.f, 5.f, 1)));

    // QObject::connect(this, &SnakeSegmentationEditor::movedOnImage, this, [&](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* event) {
    this->animate([=]() {
        snakeParameters.snake.runSegmentation(50);
        showSnakePath();
    });


    QObject::connect(this, &SnakeSegmentationEditor::clickedOnImage, this, [&](const Vector3& relPos, Vector3 value, bool leftClick, bool rightClick) {
        if (!relPos.isValid()) return;

        Vector3 pos = relPos * Vector3(100, 100, 1);
        if (snakeParameters.params->collapseFirstAndLastPoint)
            snakeParameters.snake.contour = ShapeCurve::circle(10.f, pos, 20);
        else
            snakeParameters.snake.contour = BSpline({pos - Vector3(10, 0), pos + Vector3(10, 0)}).resamplePoints(20);

        // parameters.snake.runSegmentation(10);

        showSnakePath();
    });

    QObject::connect(this, &SnakeSegmentationEditor::movedOnImage, this, [&](const Vector3& pos, const Vector3& previousPos, QMouseEvent* event) {
        bool leftPressed = event->buttons().testFlag(Qt::LeftButton);
        bool rightPressed = event->buttons().testFlag(Qt::RightButton);
        if (!leftPressed && !rightPressed) return;
        PainterToolsUI::paintImage(currentField, pos, painterParameters, rightPressed);
        currentGradientField = currentField.gradient();

        this->addImage(currentField);
    });
}

void SnakeSegmentationEditor::showSnakePath() {
    GridV3 colors(300, 300, 1);
    GridF alpha(colors.getDimensions(), 0.f);
    Vector3 ratio = (this->hasImage() ? (Vector3)colors.getDimensions() / (Vector3)this->dataModel->imageData.getImage().getDimensions() : Vector3(1.f, 1.f, 1.f));

    BSpline path;
    if (snakeParameters.params->collapseFirstAndLastPoint)
        path = ShapeCurve(snakeParameters.snake.contour).scale(ratio).getPath(100);
    else
        path = BSpline(snakeParameters.snake.contour).scale(ratio).getPath(100);


    for (size_t i = 0; i < path.size() - 1; i++) {
        PlotVectorData::drawLine(colors, Vector3::blue, path[i], path[i + 1]);
        PlotVectorData::drawLine(alpha, 1.f, path[i], path[i + 1]);
    }
    for (auto& p : snakeParameters.snake.contour) {
        colors(p * ratio) = Vector3::yellow;
        alpha(p * ratio) = 1.f;
    }

    this->setOverlay(colors, alpha, "snake");

    // std::cout << snakeParameters.snake.contour << std::endl;
    // std::cout << "Area: " << ShapeCurve(snakeParameters.snake.contour).computeArea() << "/" << snakeParameters.params->targetArea << " -- Length: " << snakeParameters.snake.contour.length() << "/" << snakeParameters.params->targetLength << std::endl;
    // this->show();
    this->draw();
}

SnakeSegmentationEditor& SnakeSegmentationEditor::associateEnvObject(EnvObject* obj)
{
    this->associatedObject = obj;
    // obj->updateFittingFunction();
    *(obj->snakeField) = SnakeImageFieldImplicit([=](const Vector3& p) { return this->currentField.interpolate(p); });
    this->snakeParameters.snake = obj->instantiate()->snake;
    return *this;
}

SnakeSegmentationEditor& SnakeSegmentationEditor::setSnakeImage(const GridF& newFieldValues)
{
    this->currentField = newFieldValues;
    this->currentGradientField = currentField.gradient();

    auto snakeField = dynamic_cast<SnakeImageFieldImplicit*>(snakeParameters.field);
    snakeField->imageField = [=](const Vector3& p) { return currentField.interpolate(p); };
    snakeField->gradientField = [=](const Vector3& p) { return currentGradientField.interpolate(p); };

    this->addImage(currentField);
    return *this;
}

SnakeSegmentationEditor& SnakeSegmentationEditor::updateToolsInterface()
{
    this->toolsInterface->clear();

    auto UI = new InterfaceUI();

    auto connectivityCostInput = new FloatInputElement("Connectivity", snakeParameters.params->connectivityCost);

    auto curvatureCostInput = new FloatInputElement("curvatureCost", snakeParameters.params->curvatureCost);
    auto imageCostInput = new FloatInputElement("imageCost", snakeParameters.params->imageCost);
    auto areaCostInput = new FloatInputElement("areaCost", snakeParameters.params->areaCost);
    auto lengthCostInput = new FloatInputElement("lengthCost", snakeParameters.params->lengthCost);
    auto slopeCostInput = new FloatInputElement("slopeCost", snakeParameters.params->slopeCost);

    auto positionCostInput = new FloatInputElement("positionCost", snakeParameters.params->positionCost);

    auto imageBordersCoefInput = new FloatInputElement("imageBordersCoef", snakeParameters.params->imageBordersCoef);
    auto imageInsideCoefInput = new FloatInputElement("imageInsideCoef", snakeParameters.params->imageInsideCoef);

    auto targetLengthInput = new FloatInputElement("targetLength", snakeParameters.params->targetLength);
    auto targetAreaInput = new FloatInputElement("targetArea", snakeParameters.params->targetArea);

    auto nbCatapillarsInput = new SliderElement("Catapillars", 0.f, 10.f, 1.0);

    auto collapseFirstAndLastPointCheckbox = new CheckboxElement("Closed shape", snakeParameters.params->collapseFirstAndLastPoint);

    auto stepSizeSlider = new SliderElement("Step size", 0.001f, 1.f, 0.001f, snakeParameters.snake.stepSize);

    nbCatapillarsInput->setValue(snakeParameters.params->nbCatapillars)
        .setOnValueChanged([=](float newVal) { snakeParameters.params->nbCatapillars = std::round(newVal);});
    /*
    UI->add(std::vector<UIElement*>{
        PainterToolsUI::createPainterToolsUI(&painterParameters),
        connectivityCostInput,
        curvatureCostInput,
        imageCostInput,
        areaCostInput,
        lengthCostInput,
        slopeCostInput,
        positionCostInput,
        imageBordersCoefInput,
        imageInsideCoefInput,
        targetLengthInput,
        targetAreaInput,
        nbCatapillarsInput,
        collapseFirstAndLastPointCheckbox,
        stepSizeSlider
    });
    */

    this->toolsInterface->add(std::move(UI));

    return *this;
}

SnakeSegmentationEditor& SnakeSegmentationEditor::updateViewOptionsInterface()
{
    this->viewOptionsInterface->clear();
    ImageViewer::updateViewOptionsInterface();
    return *this;
}
