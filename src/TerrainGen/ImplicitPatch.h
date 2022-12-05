#ifndef IMPLICITPATCH_H
#define IMPLICITPATCH_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"

class ImplicitPatch
{
public:
    ImplicitPatch();
    ImplicitPatch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, float densityValue = 1.f);
    ImplicitPatch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, std::function<float(Vector3)> compositionFunction, float densityValue = 1.f);
    ImplicitPatch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, std::function<float(Vector3)> compositionFunction, std::function<float(Vector3)> alphaDistanceFunction, float densityValue = 1.f);
    ~ImplicitPatch();

    void deleteComposables();
    ImplicitPatch clone();

    float getMaxHeight(Vector3 position);
    float evaluate(Vector3 position);
    float evaluateFromValues(float evaluationA, float evaluationB);
    /*float get(Vector3 position) {
        float composition = this->compositionFunction(position);
        return composition * this->densityValue;
    }*/
    std::map<float, float> getDensityAtPosition(Vector3 position);

    Vector3 getPositionToEvaluateComposableA(Vector3 pos);
    Vector3 getPositionToEvaluateComposableB(Vector3 pos);

    float getAlphaValue(Vector3 pos) {
        float alpha = this->alphaDistanceFunction(pos);
        return alpha;
    }

    enum CompositionFunction {
        STACK,
        BLEND,
        REPLACE,
        NONE
    };

    bool isIdentity() { return this->compositionOperation == NONE && this->getDimensions() == Vector3(); }
    bool isOperation() { return this->compositionOperation != NONE; }
    Vector3 position;
//    Vector3 boundMin;
    Vector3 dimension;
    std::function<float(Vector3)> evalFunction;
    float densityValue;
    std::function<float(Vector3)> compositionFunction;
    Matrix3<float> distanceTransform;
    std::function<float(Vector3)> alphaDistanceFunction;

    Vector3 getDimensions() const { return this->dimension/* - this->boundMin*/; }

    Matrix3<float> _cachedHeights;

    std::string toString() { return this->name + " #" + std::to_string(this->index); }

    std::string name = "Primitive";
    int index = -1;

    ImplicitPatch* composableA = nullptr;
    ImplicitPatch* composableB = nullptr;
    CompositionFunction compositionOperation = NONE;


    static ImplicitPatch createStack(ImplicitPatch *P1, ImplicitPatch *P2);
    static ImplicitPatch createReplacement(ImplicitPatch* P1, ImplicitPatch* P2);
    static ImplicitPatch createBlending(ImplicitPatch* P1, ImplicitPatch* P2);


    static int currentMaxIndex;

    static std::function<float(Vector3)> createSphereFunction(float radius);
    static std::function<float(Vector3)> createBlockFunction(float width, float depth, float height);
    static std::function<float(Vector3)> createGaussianFunction(float sigma, float width, float depth, float height);
//    static std::function<float(Vector3)> ...;

    static std::function<float(Vector3)> convert2DfunctionTo3Dfunction(std::function<float(Vector3)> func);
};
/*
class Patch2D: public ImplicitPatch
{
public:
    Patch2D();
    Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> heightFunction, float densityValue = 1.f);
    Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> heightFunction, std::function<float(Vector3)> compositionFunction, float densityValue = 1.f);
    Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> heightFunction, std::function<float(Vector3)> compositionFunction, std::function<float(Vector3)> alphaDistanceFunction, float densityValue = 1.f);
//    ~Patch2D() {}

    float getMaxHeight(Vector3 position);
    float evaluate(Vector3 position);


    Vector3 boundMin;
};

class Patch3D: public ImplicitPatch
{
public:
    Patch3D();
    Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, float densityValue = 1.f);
    Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, std::function<float(Vector3)> compositionFunction, float densityValue = 1.f);
    Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, std::function<float(Vector3)> compositionFunction, std::function<float(Vector3)> alphaDistanceFunction, float densityValue = 1.f);
    Patch3D(const Patch3D& copy);
//    ~Patch3D() {}

    float getMaxHeight(Vector3 position);
    float evaluate(Vector3 position);

    static Patch3D stack(Patch3D *P1, Patch3D *P2);
    static Patch3D replace(Patch3D* P1, Patch3D* P2);
    static Patch3D blend(Patch3D* P1, Patch3D* P2);

    Vector3 boundMin;
};
*/


#endif // IMPLICITPATCH_H
