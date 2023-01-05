#ifndef IMPLICITPATCH_H
#define IMPLICITPATCH_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "Utils/json.h"

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
    void cleanCache();

    float getMaxHeight(Vector3 position);
    float getMinHeight(Vector3 position);
    float evaluate(Vector3 position);
    float evaluateFromValues(float evaluationA, float evaluationB);
    /*float get(Vector3 position) {
        float composition = this->compositionFunction(position);
        return composition * this->densityValue;
    }*/
    std::map<float, float> getDensityAtPosition(Vector3 position);

    Vector3 getPositionToEvaluateComposableA(Vector3 pos);
    Vector3 getPositionToEvaluateComposableB(Vector3 pos);

    Vector3 getMinPosition();
    Vector3 getMaxPosition();

    void defineFunctionsBasedOnPredefinedShape();

    nlohmann::json toJson() const;
    static ImplicitPatch* fromJson(nlohmann::json json_content);


    float getAlphaValue(Vector3 pos) {
        float alpha = this->alphaDistanceFunction(pos);
        return alpha;
    }

    enum CompositionFunction {
        STACK,
        BLEND,
        REPLACE,
        NEG_STACKING,
        STACK_IN_WATER,
        NONE
    };

    enum PredefinedShapes {
        Sphere,
        Block,
        Gaussian,
        Rock,
        Mountain,
        Dune,
        Basin,
        Cave,
        Arch,
        None
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
    float blendingFactor = 1.f;
    float sigmaValue = 1.f;

    PredefinedShapes shape = None;

    Vector3 getDimensions() const { return this->dimension/* - this->boundMin*/; }

    Matrix3<float> _cachedMaxHeight;
    Matrix3<float> _cachedMinHeight;

    std::string toString();

    std::string name = "Primitive";
    int index = -1;

    ImplicitPatch* composableA = nullptr;
    ImplicitPatch* composableB = nullptr;
    CompositionFunction compositionOperation = NONE;
    std::string used_json_filename = "";


    static ImplicitPatch createStack(ImplicitPatch *P1, ImplicitPatch *P2);
    static ImplicitPatch createReplacement(ImplicitPatch* P1, ImplicitPatch* P2);
    static ImplicitPatch createBlending(ImplicitPatch* P1, ImplicitPatch* P2);
    static ImplicitPatch createNegStacking(ImplicitPatch* P1, ImplicitPatch* P2);
    static ImplicitPatch createStackInWater(ImplicitPatch* P1, ImplicitPatch* P2);
    void initCachedAttributes(ImplicitPatch* composableA = nullptr, ImplicitPatch* composableB = nullptr);


    static int currentMaxIndex;

    static std::function<float(Vector3)> createSphereFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createBlockFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createGaussianFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createRockFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createMountainFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createDuneFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createBasinFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createCaveFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createArchFunction(float sigma, float width, float depth, float height);
//    static std::function<float(Vector3)> ...;

    static std::function<float(Vector3)> convert2DfunctionTo3Dfunction(std::function<float(Vector3)> func);
};

ImplicitPatch::CompositionFunction compositionOperationFromString(std::string name);
std::string stringFromCompositionOperation(ImplicitPatch::CompositionFunction operation);
ImplicitPatch::PredefinedShapes predefinedShapeFromString(std::string name);
std::string stringFromPredefinedShape(ImplicitPatch::PredefinedShapes shape);



#endif // IMPLICITPATCH_H
