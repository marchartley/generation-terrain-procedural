#ifndef IMPLICITPATCH_H
#define IMPLICITPATCH_H

class ImplicitPatch;
class ImplicitPrimitive;
class ImplicitBinaryOperator;
class ImplicitUnaryOperator;
class ImplicitNaryOperator;
class ImplicitCSG;
class UnaryOp;

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Voxel.h"
#include "Utils/json.h"
#include "TerrainGen/TerrainModel.h"

//class TerrainMaterial;


/*
class TerrainMaterial {
public:
    float density = 1.f;
    float resistanceToWater;
    float porosity;

    Vector3 color;

    operator float() const { return this->density; }
};*/

class ImplicitPatch : public TerrainModel { // abstract
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
        FIXED_POS,
        SMOOTH_ABOVE
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
        Polygon,
        ImplicitHeightmap,
        ParametricTunnel,
        Ripple,
        None
    };

    virtual float evaluate(Vector3 pos) = 0;
//    virtual std::map<TerrainMaterial, float> getMaterials(Vector3 pos) = 0;
    virtual std::map<TerrainTypes, float> getMaterials(Vector3 pos) = 0;
    virtual float getMaxHeight(Vector3 pos);
    virtual float getMinHeight(Vector3 pos);
    virtual float getMinimalHeight(AABBox BBox);
    virtual float getMaximalHeight(AABBox BBox);
    virtual float getMinimalHeight(Vector3 minBox = Vector3::min(), Vector3 maxBox = Vector3::max());
    virtual float getMaximalHeight(Vector3 minBox = Vector3::min(), Vector3 maxBox = Vector3::max());
    std::pair<float, std::map<TerrainTypes, float> > getMaterialsAndTotalEvaluation(Vector3 pos);

    virtual AABBox getSupportBBox() = 0;
    virtual AABBox getBBox() = 0;
    Vector3 getDimensions();
    Vector3 getSupportDimensions();

    Vector3 getNormal(Vector3 pos);

    void setIndex(int newIndex = -1);

    virtual void update() = 0;

    virtual std::string toString() = 0;

    virtual nlohmann::json toJson() = 0;
    static ImplicitPatch* fromJson(nlohmann::json content);

    virtual void updateCache();

    virtual ImplicitPatch* copy() const = 0;

    bool checkIsInGround(Vector3 position);

    virtual void initMap() { }

    virtual bool undo() { return false; }
    virtual bool redo() { return false; }

    virtual void saveMap(std::string filename) { }
    virtual void retrieveMap(std::string filename) { this->fromJson(nlohmann::json::parse(std::ifstream(filename))); }
    virtual Mesh getGeometry();

    virtual Vector3 getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos = Vector3(false), Vector3 maxPos = Vector3(false));

    virtual std::string toShortString() { return ""; };

    virtual float getHeight(float x, float y);
    virtual float getHeight(Vector3 pos);

    virtual bool contains(Vector3 v);
    virtual bool contains(float x, float y, float z);

    virtual bool contains(PredefinedShapes shape) = 0;
    virtual std::vector<ImplicitPatch*> findAll(PredefinedShapes shape) = 0;

    Vector3 getGlobalPositionOf(Vector3 posInsidePatch);
    Vector3 getLocalPositionOf(Vector3 globalPosition);

//    virtual size_t getCurrentHistoryIndex() const;

    virtual float getSizeX() { return this->getBBox().max().x; }
    virtual float getSizeY() { return this->getBBox().max().y; }
    virtual float getSizeZ() { return this->getBBox().max().z; }

    Matrix3<float> getVoxelized(Vector3 dimensions = Vector3(false), Vector3 scale = Vector3(1.f, 1.f, 1.f));

    int index = -1;
    std::string name;
    bool mirrored = false;

    std::string used_json_filename = "";

    BSpline optionalCurve;

    static ImplicitPatch* createIdentity();
    static ImplicitPrimitive *createPredefinedShape(PredefinedShapes shape, Vector3 dimensions, float additionalParam, BSpline parametricCurve = BSpline());
    static std::function<float(Vector3)> createPredefinedShapeFunction(PredefinedShapes shape, Vector3 dimensions, float additionalParam, BSpline parametricCurve = BSpline());
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
    static std::function<float(Vector3)> createPolygonFunction(float sigma, float width, float depth, float height, BSpline path);
    static std::function<float(Vector3)> createParametricTunnelFunction(float sigma, float width, float depth, float height, BSpline path);
    static std::function<float(Vector3)> createRippleFunction(float sigma, float width, float depth, float height);
    static std::function<float(Vector3)> createIdentityFunction(float sigma, float width, float depth, float height);
//    static std::function<float(Vector3)> ...;

    static std::function<float(Vector3)> convert2DfunctionTo3Dfunction(std::function<float(Vector3)> func);

    ImplicitNaryOperator *getParent() const;
    void setParent(ImplicitNaryOperator* newParent);
    ImplicitNaryOperator* parent = nullptr;

    static int currentMaxIndex;

    static std::string json_identifier;

//protected:
    Matrix3<float> _cachedMinHeight;
    Matrix3<float> _cachedMaxHeight;
    Matrix3<float> _cachedVoxelized;
    bool _cached = false;
};

class ImplicitPrimitive : public ImplicitPatch {
public:
    ImplicitPrimitive();

    float evaluate(Vector3 pos);
    std::map<TerrainTypes, float> getMaterials(Vector3 pos);

    AABBox getSupportBBox();
    AABBox getBBox();
    void update();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    void setDimensions(Vector3 newDimensions);
    void setSupportDimensions(Vector3 newSupportDimensions);

    bool contains(PredefinedShapes shape) { return this->predefinedShape == shape; }
    std::vector<ImplicitPatch*> findAll(PredefinedShapes shape);

    virtual ImplicitPatch* copy() const;

    Vector3 position = Vector3(false);
    Vector3 dimensions = Vector3(false);
    Vector3 supportDimensions = Vector3(false);
    std::function<float(Vector3)> evalFunction;
    TerrainTypes material = WATER;

    PredefinedShapes predefinedShape = None;
    std::vector<float> parametersProvided;

    std::string heightmapFilename = "";
    Matrix3<float> cachedHeightmap;

    static ImplicitPrimitive* fromHeightmap(std::string filename, Vector3 dimensions = Vector3(false), ImplicitPrimitive *prim = nullptr);
    static ImplicitPrimitive* fromHeightmap(Matrix3<float> heightmap, std::string filename = "", ImplicitPrimitive *prim = nullptr);
};

class ImplicitNaryOperator : public ImplicitPatch {
public:
    ImplicitNaryOperator();

    float evaluate(Vector3 pos);
    std::map<TerrainTypes, float> getMaterials(Vector3 pos);
    AABBox getSupportBBox();
    AABBox getBBox();
    void update();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    bool contains(PredefinedShapes shape);
    virtual std::vector<ImplicitPatch*> findAll(PredefinedShapes shape);

    void augment();

    void updateCache();

    virtual ImplicitPatch* copy() const;

    virtual void addChild(ImplicitPatch* newChild, int index = -1);
    virtual void deleteAllChildren();

    std::vector<ImplicitPatch*> composables;
};

class ImplicitBinaryOperator : public ImplicitNaryOperator {
public:
    ImplicitBinaryOperator();

    float evaluate(Vector3 pos);
    std::map<TerrainTypes, float> getMaterials(Vector3 pos);
    float evaluateFromAandB(float evalA, float evalB);

    float evaluateA(Vector3 pos);
    float evaluateB(Vector3 pos);

    std::map<TerrainTypes, float> getMaterialsA(Vector3 pos);
    std::map<TerrainTypes, float> getMaterialsB(Vector3 pos);
    std::pair<float, std::map<TerrainTypes, float>> getMaterialsAndTotalEvaluationA(Vector3 pos);
    std::pair<float, std::map<TerrainTypes, float>> getMaterialsAndTotalEvaluationB(Vector3 pos);

    bool contains(PredefinedShapes shape);
    virtual std::vector<ImplicitPatch*> findAll(PredefinedShapes shape);

    AABBox getSupportBBox();
    AABBox getBBox();
    void update();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    void updateCache();

    void swapAB();


    Vector3 getEvaluationPositionForComposableA(Vector3 pos);
    Vector3 getEvaluationPositionForComposableB(Vector3 pos);

    virtual ImplicitPatch* copy() const;
    virtual void addChild(ImplicitPatch* newChild, int index);
    virtual void deleteAllChildren();

//    ImplicitPatch* composableA = nullptr;
//    ImplicitPatch* composableB = nullptr;
    ImplicitPatch*& composableA();
    ImplicitPatch*& composableB();

    CompositionFunction composeFunction;
    PositionalLabel positionalB;

    float blendingFactor = 2.f; // To replace with a full function (or curve)

    bool withIntersectionOnB = false; // Intersection or Union
};


class ImplicitUnaryOperator : public ImplicitNaryOperator {
public:
    using TransformFunction = std::function<Vector3(const Vector3&)>;

    ImplicitUnaryOperator();

    void addTransformation(TransformFunction transform, TransformFunction inverse) {
        transformations.emplace_back(std::move(transform), std::move(inverse));
    }

    float evaluate(Vector3 point);

    void addTranslation(const Vector3& translation);

    void addRotation(Vector3 rotationAngles);

    void addScaling(const Vector3& scaleFactors);

    void addDisplacementField(const Matrix3<Vector3>& displacementField);

    Vector3 applyTransform(Vector3 pos) const;
    Vector3 inverseTransform(Vector3 pos) const;


    std::map<TerrainTypes, float> getMaterials(Vector3 pos);

    AABBox getSupportBBox();
    AABBox getBBox();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    bool contains(PredefinedShapes shape);
    virtual std::vector<ImplicitPatch*> findAll(PredefinedShapes shape);

//    std::function<Vector3(Vector3)> wrapFunction;
//    std::function<Vector3(Vector3)> unwrapFunction;
    std::function<float(Vector3)> noiseFunction;
    std::vector<UnaryOp> transforms;

//    virtual ImplicitPatch* copy() const;
    virtual void addChild(ImplicitPatch* newChild, int index = 0);
    virtual void deleteAllChildren();

    ImplicitPatch*& composableA();

    void translate(Vector3 translation);
    void rotate(Vector3 angles); //float angleX, float angleY, float angleZ);
    void scale(Vector3 scaleFactor);

    void addRandomNoise(float amplitude, float period = 20.f, float offset = 10.f);
    void addRandomWrap(float amplitude, float period = 20.f, float offset = 10.f);
    void addWrapFunction(Matrix3<Vector3> func);
    void spread(float factor = 1.f);
    void addWavelets();

    Vector3 _translation = Vector3(0, 0, 0); // Should not be here, just used to be stored in files
    Vector3 _rotation = Vector3(0, 0, 0); // Should not be here, just used to be stored in files
    Vector3 _scale = Vector3(1, 1, 1);
    Vector3 _distortion = Vector3(0.f, 0.f, 0.f);
    Vector3 _noise = Vector3(0.f, 0.f, 0.f);
    float _spreadingFactor = 0.f;

private:
    std::vector<std::pair<TransformFunction, TransformFunction>> transformations;
};

class ImplicitTranslation : public ImplicitUnaryOperator {
public:
    ImplicitTranslation() : ImplicitUnaryOperator() { this->name = "Translation"; }

    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);
};

class ImplicitRotation : public ImplicitUnaryOperator {
public:
    ImplicitRotation() : ImplicitUnaryOperator() { this->name = "Rotation"; }

    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);
};

class ImplicitScaling : public ImplicitUnaryOperator {
public:
    ImplicitScaling() : ImplicitUnaryOperator() { this->name = "Scaling"; }

    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);
};

class ImplicitWraping : public ImplicitUnaryOperator {
public:
    ImplicitWraping() : ImplicitUnaryOperator() { this->name = "Wraping"; }

    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);
};

class ImplicitNoise : public ImplicitUnaryOperator {
public:
    ImplicitNoise() : ImplicitUnaryOperator() { this->name = "Noise"; }

    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);
};

class ImplicitSpread : public ImplicitUnaryOperator {
public:
    ImplicitSpread() : ImplicitUnaryOperator() { this->name = "Spread"; }

    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);
};




/*
class ImplicitUnaryOperator : public ImplicitNaryOperator {
public:
    ImplicitUnaryOperator();

    float evaluate(Vector3 pos);
    std::map<TerrainTypes, float> getMaterials(Vector3 pos);

    AABBox getSupportBBox();
    AABBox getBBox();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    bool contains(PredefinedShapes shape);
    virtual std::vector<ImplicitPatch*> findAll(PredefinedShapes shape);

    std::function<Vector3(Vector3)> wrapFunction;
    std::function<Vector3(Vector3)> unwrapFunction;
    std::function<float(Vector3)> noiseFunction;
    std::vector<UnaryOp> transforms;

//    virtual ImplicitPatch* copy() const;
    virtual void addChild(ImplicitPatch* newChild, int index = 0);
    virtual void deleteAllChildren();

    ImplicitPatch*& composableA();

    void translate(Vector3 translation);
    void rotate(float angleX, float angleY, float angleZ);
    void scale(Vector3 scaleFactor);

    void addRandomNoise(float amplitude, float period = 20.f, float offset = 10.f);
    void addRandomWrap(float amplitude, float period = 20.f, float offset = 10.f);
    void addWrapFunction(Matrix3<Vector3> func);
    void spread(float factor = 1.f);
    void addWavelets();

    Vector3 _translation = Vector3(0, 0, 0); // Should not be here, just used to be stored in files
    Vector3 _rotation = Vector3(0, 0, 0); // Should not be here, just used to be stored in files
    Vector3 _scale = Vector3(1, 1, 1);
    Vector3 _distortion = Vector3(0.f, 0.f, 0.f);
    Vector3 _noise = Vector3(0.f, 0.f, 0.f);
    float _spreadingFactor = 0.f;
};
*/







class UnaryOp {
public:
    UnaryOp();
    std::function<Vector3(Vector3)> wrap;
    std::function<Vector3(Vector3)> unwrap;
};
class UnaryOpTranslate: public UnaryOp {
public:
    UnaryOpTranslate(Vector3 translation);
};
class UnaryOpRotate: public UnaryOp {
public:
    UnaryOpRotate(Vector3 rotationAngles, Vector3 center);
};
class UnaryOpScale: public UnaryOp {
public:
    UnaryOpScale(Vector3 scaling, Vector3 center);
};
class UnaryOpWrap: public UnaryOp {
public:
    UnaryOpWrap(FastNoiseLite noise, Vector3 strength);
    UnaryOpWrap(std::function<Vector3(Vector3)> func);
    UnaryOpWrap(Matrix3<Vector3> wrapper);
    Matrix3<Vector3> wrapper;
};
class UnaryOpSpread: public UnaryOp {
public:
    UnaryOpSpread(AABBox BBox, float spreadFactor);
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

ImplicitPatch::CompositionFunction compositionOperationFromString(std::string name);
std::string stringFromCompositionOperation(ImplicitPatch::CompositionFunction operation);
ImplicitPatch::PositionalLabel positionalLabelFromString(std::string name);
std::string stringFromPositionalLabel(ImplicitPatch::PositionalLabel label);
ImplicitPatch::PredefinedShapes predefinedShapeFromString(std::string name);
std::string stringFromPredefinedShape(ImplicitPatch::PredefinedShapes shape);
TerrainTypes materialFromString(std::string name);
std::string stringFromMaterial(TerrainTypes material);


#endif // IMPLICITPATCH_H
