#ifndef KARSTPATHGENERATIONINTERFACE_H
#define KARSTPATHGENERATIONINTERFACE_H

#include "Interface/ControlPoint.h"
#include "Interface/InteractiveVector.h"
#include "Interface/Slider3D.h"
#include "Karst/KarstPathsGeneration.h"
#include "Graphics/Mesh.h"

#include <QWidget>

class KarstPathGenerationInterface : public QWidget
{
    Q_OBJECT
public:
    KarstPathGenerationInterface();
    KarstPathGenerationInterface(KarstPathsGeneration* karstCreator, Vector3 AABBoxMinPos, Vector3 AABBoxMaxPos);

    void display();

    InteractiveVector *fractureVector;
    Slider3D *waterHeightSlider;
    Mesh waterLevelMesh;

public Q_SLOTS:
    void updateFracture(Vector3 newFractureDir = Vector3());
    void updateWaterHeight(float newHeight = 0.f);

public:
    Vector3 AABBoxMinPos;
    Vector3 AABBoxMaxPos;

    KarstPathsGeneration* karstCreator = nullptr;
};

#endif // KARSTPATHGENERATIONINTERFACE_H
