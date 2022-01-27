#include "KarstPathGenerationInterface.h"

KarstPathGenerationInterface::KarstPathGenerationInterface()
    : karstCreator(nullptr)
{
}

KarstPathGenerationInterface::KarstPathGenerationInterface(KarstPathsGeneration* karstCreator, Vector3 AABBoxMinPos, Vector3 AABBoxMaxPos)
    : AABBoxMinPos(AABBoxMinPos), AABBoxMaxPos(AABBoxMaxPos)
{
    this->karstCreator = karstCreator;

    this->fractureVector = new InteractiveVector(AABBoxMinPos, Vector3(AABBoxMinPos.x, AABBoxMinPos.y, AABBoxMaxPos.z));
    this->waterHeightSlider = new Slider3D(Vector3(AABBoxMinPos.x - 10, AABBoxMinPos.y - 10, AABBoxMinPos.z), Vector3(AABBoxMinPos.x - 10, AABBoxMinPos.y - 10, AABBoxMaxPos.z), 0.f, AABBoxMinPos.z, AABBoxMaxPos.z);
//    this->fractureVector.onUpdate([&]() -> void {this->updateFracture(); });
//    this->fractureVector.onUpdate(std::bind(&KarstPathGenerationInterface::updateFracture, this));
//    this->updateFracture();

    QObject::connect(this->fractureVector, &InteractiveVector::modified, this, &KarstPathGenerationInterface::updateFracture);
    QObject::connect(this->waterHeightSlider, &Slider3D::valueChanged, this, &KarstPathGenerationInterface::updateWaterHeight);

    Vector3 offX = Vector3(AABBoxMaxPos.x - AABBoxMinPos.x, 0, 0);
    Vector3 offY = Vector3(0, AABBoxMaxPos.y - AABBoxMinPos.y, 0);
    Vector3 offZ = Vector3(0, 0, this->waterHeightSlider->sliderControlPoint->position.z);
    this->waterLevelMesh.fromArray({
                                       {AABBoxMinPos + offZ, AABBoxMinPos + offX + offZ,
                                        AABBoxMinPos + offX + offZ, AABBoxMinPos + offX + offY + offZ,
                                        AABBoxMinPos + offX + offY + offZ, AABBoxMinPos + offY + offZ,
                                        AABBoxMinPos + offY + offZ, AABBoxMinPos + offZ}
                                   });
}

void KarstPathGenerationInterface::display()
{
//    this->updateFracture(); // Should be in the callbacks, but I didn't manage to do it...
//    this->updateWaterHeight();
    this->fractureVector->display();
    this->waterLevelMesh.display(GL_LINES);
    this->waterHeightSlider->display();
}

void KarstPathGenerationInterface::updateFracture(Vector3 newFractureDir)
{
    if (karstCreator->fracturesDirections.empty())
        this->karstCreator->addFractureDirection(FractureDirection(Vector3(0, 0, 0), 1.f));
    karstCreator->fracturesDirections.back().direction = newFractureDir;
}

void KarstPathGenerationInterface::updateWaterHeight(float newHeight)
{
    if (karstCreator->waterHeights.empty())
        this->karstCreator->addWaterHeight(WaterHeight(0.f, 1.f));
    karstCreator->waterHeights.back().height = newHeight;
    Vector3 offX = Vector3(AABBoxMaxPos.x - AABBoxMinPos.x, 0, 0);
    Vector3 offY = Vector3(0, AABBoxMaxPos.y - AABBoxMinPos.y, 0);
    Vector3 offZ = Vector3(0, 0, this->waterHeightSlider->sliderControlPoint->position.z);
    this->waterLevelMesh.fromArray({
                                       {AABBoxMinPos + offZ, AABBoxMinPos + offX + offZ,
                                        AABBoxMinPos + offX + offZ, AABBoxMinPos + offX + offY + offZ,
                                        AABBoxMinPos + offX + offY + offZ, AABBoxMinPos + offY + offZ,
                                        AABBoxMinPos + offY + offZ, AABBoxMinPos + offZ}
                                   });

}
