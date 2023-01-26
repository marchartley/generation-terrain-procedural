#ifndef IMPLICITPATCH_H
#define IMPLICITPATCH_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Voxel.h"
#include "Utils/json.h"

//#define useBefore

#ifndef useBefore

//class TerrainMaterial;

class ImplicitPatch;
class ImplicitPrimitive;
class ImplicitOperator;
class ImplicitUnaryOperator;
class ImplicitCSG;

/*
class TerrainMaterial {
public:
    float density = 1.f;
    float resistanceToWater;
    float porosity;

    Vector3 color;

    operator float() const { return this->density; }
};*/

class ImplicitPatch { // abstract
public:
    ImplicitPatch();
    virtual ~ImplicitPatch() = default;

    enum CompositionFunction {
        STACK,
        BLEND,
        REPLACE,
        ONE_SIDE_BLEND,
        NONE
    };

    enum PositionalLabel {
        ABOVE,
        INSIDE_BOTTOM,
        INSIDE_TOP,
        FIXED_POS
    };

    enum PredefinedShapes {
        Sphere,
        Block,
        Gaussian,
        Cylinder,
        Rock,
        Mountain,
        Dune,
        Basin,
        Cave,
        Arch,
        Noise2D,
        MountainChain,
        None
    };

    virtual float evaluate(Vector3 pos) = 0;
//    virtual std::map<TerrainMaterial, float> getMaterials(Vector3 pos) = 0;
    virtual std::map<TerrainTypes, float> getMaterials(Vector3 pos) = 0;
    virtual float getMaxHeight(Vector3 pos);
    virtual float getMinHeight(Vector3 pos);
    std::pair<float, std::map<TerrainTypes, float> > getMaterialsAndTotalEvaluation(Vector3 pos);

    virtual std::pair<Vector3, Vector3> getSupportBBox() = 0;
    virtual std::pair<Vector3, Vector3> getBBox() = 0;
    Vector3 getDimensions();
    Vector3 getSupportDimensions();

    Vector3 getNormal(Vector3 pos);

    void setIndex(int newIndex = -1);

    virtual void update() = 0;

    virtual std::string toString() = 0;

    virtual nlohmann::json toJson() = 0;
    static ImplicitPatch* fromJson(nlohmann::json content);

    void updateCache();

    int index = -1;
    std::string name;

    std::string used_json_filename = "";

    BSpline optionalCurve;

    static ImplicitPatch* createPredefinedShape(PredefinedShapes shape, Vector3 dimensions, float additionalParam, BSpline parametricCurve = BSpline());
    static float isovalue;
    static float zResolution;


    static std::function<float(Vector3)> createSphereFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createBlockFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createGaussianFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createCylinderFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createRockFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createMountainFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createDuneFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createBasinFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createCaveFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createArchFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createNoise2DFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createMountainChainFunction(float sigma, float width, float depth, float height, BSpline path);
    static std::function<float(Vector3)> createIdentityFunction(float sigma, float width, float depth, float height);
//    static std::function<float(Vector3)> ...;

    static std::function<float(Vector3)> convert2DfunctionTo3Dfunction(std::function<float(Vector3)> func);

    static int currentMaxIndex;

    static std::string json_identifier;

protected:
    Matrix3<float> _cachedMinHeight;
    Matrix3<float> _cachedMaxHeight;
};

class ImplicitPrimitive : public ImplicitPatch {
public:
    ImplicitPrimitive();

    float evaluate(Vector3 pos);
    std::map<TerrainTypes, float> getMaterials(Vector3 pos);

    std::pair<Vector3, Vector3> getSupportBBox();
    std::pair<Vector3, Vector3> getBBox();
    void update();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    void setDimensions(Vector3 newDimensions);
    void setSupportDimensions(Vector3 newSupportDimensions);

    Vector3 position = Vector3(false);
    Vector3 dimensions = Vector3(false);
    Vector3 supportDimensions = Vector3(false);
    std::function<float(Vector3)> evalFunction;
    TerrainTypes material = WATER;

    PredefinedShapes predefinedShape = None;
    std::vector<float> parametersProvided;
};

class ImplicitOperator : public ImplicitPatch {
public:
    ImplicitOperator();

    float evaluate(Vector3 pos);
    std::map<TerrainTypes, float> getMaterials(Vector3 pos);
    float evaluateFromAandB(float evalA, float evalB);

    float evaluateA(Vector3 pos);
    float evaluateB(Vector3 pos);

    std::map<TerrainTypes, float> getMaterialsA(Vector3 pos);
    std::map<TerrainTypes, float> getMaterialsB(Vector3 pos);
    std::pair<float, std::map<TerrainTypes, float>> getMaterialsAndTotalEvaluationA(Vector3 pos);
    std::pair<float, std::map<TerrainTypes, float>> getMaterialsAndTotalEvaluationB(Vector3 pos);

    std::pair<Vector3, Vector3> getSupportBBox();
    std::pair<Vector3, Vector3> getBBox();
    void update();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    void updateCache();

    Vector3 getEvaluationPositionForComposableA(Vector3 pos);
    Vector3 getEvaluationPositionForComposableB(Vector3 pos);

    ImplicitPatch* composableA = nullptr;
    ImplicitPatch* composableB = nullptr;

    CompositionFunction composeFunction;
    PositionalLabel positionalB;

    float blendingFactor = 2.f; // To replace with a full function (or curve)

    bool withIntersectionOnB = false; // Intersection or Union
};

class ImplicitUnaryOperator : public ImplicitOperator {
public:
    ImplicitUnaryOperator();

    float evaluate(Vector3 pos);
    std::map<TerrainTypes, float> getMaterials(Vector3 pos);

    std::pair<Vector3, Vector3> getSupportBBox();
    std::pair<Vector3, Vector3> getBBox();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    std::function<Vector3(Vector3)> wrapFunction;
    std::function<float(Vector3)> noiseFunction;

    void translate(Vector3 translation);
    void rotate(float angleX, float angleY, float angleZ);
    void scale(Vector3 scaleFactor);

    void addRandomNoise(float amplitude, float period = 20.f, float offset = 10.f);
    void addRandomWrap(float amplitude, float period = 20.f, float offset = 10.f);

    Vector3 _translation = Vector3(0, 0, 0); // Should not be here, just used to be stored in files
    Vector3 _rotation = Vector3(0, 0, 0); // Should not be here, just used to be stored in files
    Vector3 _scale = Vector3(1, 1, 1);
    Vector3 _distortion = Vector3(0.f, 0.f, 0.f);
    Vector3 _noise = Vector3(0.f, 0.f, 0.f);
};

/*
class ImplicitCSG : public ImplicitOperator {
public:
    ImplicitCSG();

    enum CSG_Operation {
        Union,
        Difference,
        Intersection
    };

    float evaluate(Vector3 pos);
    std::map<TerrainMaterial, float> getMaterials(Vector3 pos);
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    CSG_Operation operation = Union;
};
*/

#else
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
        std::map<float, float> getDensityAtPosition(Vector3 position, Vector3 worldPosition = Vector3(false));

        Vector3 getPositionToEvaluateComposableA(Vector3 pos);
        Vector3 getPositionToEvaluateComposableB(Vector3 pos);

        Vector3 getMinPosition();
        Vector3 getMaxPosition();

        void defineFunctionsBasedOnPredefinedShape();

        nlohmann::json toJson();
        static ImplicitPatch* fromJson(nlohmann::json json_content);


        float getAlphaValue(Vector3 pos) {
            float alpha = this->alphaDistanceFunction(pos);
            return alpha;
        }

        enum CompositionFunction {
            STACK,
            BLEND,
            REPLACE,
    //        NEG_STACKING,
    //        STACK_IN_WATER,
            NONE
        };

        enum PositionalLabel {
            ABOVE,
            INSIDE_BOTTOM,
            INSIDE_TOP,
            FIXED
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

        Vector3 getDimensions() { return /*this->getMaxPosition() - this->getMinPosition()*/this->dimension/* - this->boundMin*/; }

        Matrix3<float> _cachedMaxHeight;
        Matrix3<float> _cachedMinHeight;

        std::string toString();

        std::string name = "Primitive";
        int index = -1;

        ImplicitPatch* composableA = nullptr;
        ImplicitPatch* composableB = nullptr;
        CompositionFunction compositionOperation = NONE;
        PositionalLabel positioning = FIXED;
        std::string used_json_filename = "";


        static ImplicitPatch createStack(ImplicitPatch *P1, ImplicitPatch *P2, PositionalLabel pos);
        static ImplicitPatch createReplacement(ImplicitPatch* P1, ImplicitPatch* P2, PositionalLabel pos);
        static ImplicitPatch createBlending(ImplicitPatch* P1, ImplicitPatch* P2, PositionalLabel pos);
    //    static ImplicitPatch createNegStacking(ImplicitPatch* P1, ImplicitPatch* P2);
    //    static ImplicitPatch createStackInWater(ImplicitPatch* P1, ImplicitPatch* P2);
        void initCachedAttributes(ImplicitPatch* composableA = nullptr, ImplicitPatch* composableB = nullptr, PositionalLabel pos = FIXED);


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

#endif


ImplicitPatch::CompositionFunction compositionOperationFromString(std::string name);
std::string stringFromCompositionOperation(ImplicitPatch::CompositionFunction operation);
ImplicitPatch::PositionalLabel positionalLabelFromString(std::string name);
std::string stringFromPositionalLabel(ImplicitPatch::PositionalLabel label);
ImplicitPatch::PredefinedShapes predefinedShapeFromString(std::string name);
std::string stringFromPredefinedShape(ImplicitPatch::PredefinedShapes shape);
TerrainTypes materialFromString(std::string name);
std::string stringFromMaterial(TerrainTypes material);


#endif // IMPLICITPATCH_H
