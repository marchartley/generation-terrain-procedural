#ifndef INTERACTIVEVECTOR_H
#define INTERACTIVEVECTOR_H

#include "DataStructure/Vector3.h"
#include "Interface/ControlPoint.h"
#include "Graphics/Mesh.h"

class InteractiveVector : public QObject
{
    Q_OBJECT
public:
    InteractiveVector();
    InteractiveVector(Vector3 end);
    InteractiveVector(Vector3 start, Vector3 end);

    void display();

Q_SIGNALS:
    void modified(Vector3 newVector);
    void startingModified(Vector3 newVector);
    void endingModified(Vector3 newVector);

public:
    Vector3 getResultingVector() { return getEndingVector() - getStartingVector(); }
    Vector3 getStartingVector() { return this->startingControlPoint->position; }
    Vector3 getEndingVector() { return this->endingControlPoint->position; }

    ControlPoint *startingControlPoint;
    ControlPoint *endingControlPoint;
    Mesh arrowMesh;

    std::vector<Vector3> getArrowPath();

    static std::shared_ptr<Shader> base_shader;
};

#endif // INTERACTIVEVECTOR_H
