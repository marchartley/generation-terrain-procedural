#ifndef KARSTPATHGENERATIONINTERFACE_H
#define KARSTPATHGENERATIONINTERFACE_H

#include "Interface/ControlPoint.h"
#include "Interface/InteractiveVector.h"
#include "Interface/Slider3D.h"
#include "Karts/KarstPathsGeneration.h"
#include "Graphics/Mesh.h"

class KarstPathGenerationInterface
{
public:
    KarstPathGenerationInterface();
    KarstPathGenerationInterface(KarstPathsGeneration* karstCreator, Vector3 AABBoxMinPos, Vector3 AABBoxMaxPos);

    void display();

    InteractiveVector fractureVector;
    Slider3D waterHeightSlider;
    Mesh waterLevelMesh;

    void updateFracture();
    void updateWaterHeight();

    Vector3 AABBoxMinPos;
    Vector3 AABBoxMaxPos;

    KarstPathsGeneration* karstCreator = nullptr;
};

#endif // KARSTPATHGENERATIONINTERFACE_H
