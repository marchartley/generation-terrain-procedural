#include "GanUIPainter.h"

GanUIPainter::GanUIPainter(const std::string& name, QWidget *parent) : AbstractPlotter(name, name, parent)
{
    auto brushSizeSlider = new SliderElement("Brush size", 1, 50, 1, this->brushSize);
    // auto colorSlider = std::make_shared<SliderElement>("Biome", .2f, 1.f, 0.01f, this->colorIndex); // Ignore the 0-0.2 range as it is a transition from abyss to island center.
    auto colorSlider = new SliderElement("Biome", 0, 4, 1, this->colorIndex);
    // auto sharpnessSlider = std::make_shared<SliderElement>("Sharpness", 0.01f, 10.f, 0.01f, this->sharpness);
    colorSlider->slider()->addTicks({
        {0, "Island"},
        {1, "Beach"},
        {2, "Lagoon"},
        {3, "Reef"},
        {4, "Abyss"}
    });

    this->toolsInterface->add({
        brushSizeSlider,
        colorSlider,
        // sharpnessSlider
    });

    this->addImage(GridV3(256, 256, 1, Vector3(1.f, 0.f, 0.f)));
    QObject::connect(this, &GanUIPainter::clickedOnImage, this, [&](const Vector3& pos, Vector3 value) {
        if (pos.isValid()) {
            Vector3i p = pos * this->dataModel->getImage().getDimensions();
            this->drawStroke(p, p);
            this->draw();
        }
    });
    QObject::connect(this, &GanUIPainter::movedOnImage, this, [&](const Vector3& pos, const Vector3& previousPos, QMouseEvent* event) {
        this->currentColor = this->getColorFromIndex(this->colorIndex);
        if (event->buttons().testFlag(Qt::LeftButton)) {
            displayProcessTime("Brushing... ", [&]() {
                this->drawStroke(previousPos, pos);
                this->draw();
            }, false);
        }
    });
}

void GanUIPainter::drawStroke(const Vector3& pStart, const Vector3& pEnd)
{
    auto currentImage = this->dataModel->getImage();
    Vector3i brushCenter = Vector3i(brushSize / 2.f, brushSize / 2.f);
    Vector3 minMask = Vector3::min(pStart, pEnd) - brushCenter;
    Vector3 maxMask = Vector3::max(pStart, pEnd) + brushCenter;
    Vector3 modifsDimensions = maxMask - minMask;
    GridF mask = GridF(modifsDimensions.x(), modifsDimensions.y(), 1);
    Vector3 p1 = pStart - minMask;
    Vector3 p2 = pEnd - minMask;
    mask.iterateParallel([&](const Vector3i& p) {
        mask[p] = std::min(1.f, (p - Collision::projectPointOnSegment(p, p1, p2)).length() / (brushSize * .5f));
    });
    mask.iterateParallel([&](const Vector3i& p) {
        currentImage[p + minMask] = (mask[p] >= 1.f ? currentImage[p + minMask] : currentColor); //Vector3::slerp(interpolation::power_wyvill(mask[p], std::max(0.001f, this->sharpness)), currentImage[p + minMask], currentColor);
    });
    this->dataModel->addImage(currentImage);
}

GanUIPainter& GanUIPainter::updateUI()
{
    AbstractPlotter::updateUI();
    return *this;
}

Vector3 GanUIPainter::getColorFromIndex(float colorIndex) {
    // return HSVtoRGB(this->colorIndex, 1.f, 1.f);
    // return colorPalette(this->colorIndex, {Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1)});
    // return Vector3::slerp(this->colorIndex, Vector3(1, 0, 0), Vector3(0, 0, 1));
    switch(int(colorIndex)) {
    case 0: return Vector3(1, 1, 0); // Island = yellow
    case 1: return Vector3(0, 1, 0); // Beach = green
    case 2: return Vector3(0, 1, 1); // Lagoon = cyan
    case 3: return Vector3(0, 0, 1); // Reef = blue
    case 4: return Vector3(1, 0, 0); // Abyss = red
    default: return Vector3(0, 0, 0);
    }
}
