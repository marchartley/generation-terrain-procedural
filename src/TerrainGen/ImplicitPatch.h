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
        DistanceMap,
        None
    };

    virtual float evaluate(const Vector3& pos) const = 0;
//    virtual std::map<TerrainMaterial, float> getMaterials(const Vector3& pos) = 0;
    virtual std::map<TerrainTypes, float> getMaterials(const Vector3& pos) const = 0;
    virtual float getMaxHeight(const Vector3& pos);
    virtual float getMinHeight(const Vector3& pos);
    virtual float getMinimalHeight(AABBox BBox);
    virtual float getMaximalHeight(AABBox BBox);
    virtual float getMinimalHeight(const Vector3 &minBox = Vector3::min(), const Vector3& maxBox = Vector3::max());
    virtual float getMaximalHeight(const Vector3& minBox = Vector3::min(), const Vector3& maxBox = Vector3::max());
    std::pair<float, std::map<TerrainTypes, float> > getMaterialsAndTotalEvaluation(const Vector3& pos) const;
    std::pair<float, std::map<TerrainTypes, float> > getBinaryMaterialsAndTotalEvaluation(const Vector3& pos) const;

    virtual AABBox getSupportBBox() const = 0;
    virtual AABBox getBBox() const = 0;
    Vector3 getDimensions() const;
    Vector3 getSupportDimensions() const;

    virtual GridV3 getNormals();
    virtual Vector3 getNormal(const Vector3& pos) const;

    void setIndex(int newIndex = -1);

    virtual void update() = 0;

    virtual std::string toString() = 0;

    virtual nlohmann::json toJson() = 0;
    static ImplicitPatch* fromJson(nlohmann::json content);

    virtual void updateCache();

    virtual ImplicitPatch* copy() const = 0;

    bool checkIsInGround(const Vector3& position);

    virtual void initMap() { }

    virtual bool undo() { return false; }
    virtual bool redo() { return false; }

    virtual void saveMap(std::string filename) { }
    virtual void retrieveMap(std::string filename) { this->fromJson(nlohmann::json::parse(std::ifstream(filename))); }
    virtual Mesh getGeometry(const Vector3& dimensions = Vector3(false));

    virtual Vector3 getIntersection(const Vector3& origin, const Vector3& dir, const Vector3 &minPos = Vector3(false), const Vector3 &maxPos = Vector3(false));

    virtual std::string toShortString() { return ""; }

    virtual float getHeight(float x, float y);
    virtual float getHeight(const Vector3& pos);

    virtual bool contains(const Vector3& v);
    virtual bool contains(float x, float y, float z);

    virtual bool contains(PredefinedShapes shape) = 0;
    virtual std::vector<ImplicitPatch*> findAll(PredefinedShapes shape) = 0;

    Vector3 getGlobalPositionOf(const Vector3& posInsidePatch);
    Vector3 getLocalPositionOf(const Vector3& globalPosition);

//    virtual size_t getCurrentHistoryIndex() const;

    virtual float getSizeX() const { return this->getBBox().max().x; }
    virtual float getSizeY() const { return this->getBBox().max().y; }
    virtual float getSizeZ() const { return this->getBBox().max().z; }

    GridF getVoxelized(const Vector3& dimensions = Vector3(false), const Vector3& scale = Vector3(1.f, 1.f, 1.f));

    int index = -1;
    std::string name;
    bool mirrored = false;

    std::string used_json_filename = "";

    BSpline optionalCurve;

    static ImplicitPatch* createIdentity();
    static ImplicitPrimitive *createPredefinedShape(PredefinedShapes shape, const Vector3& dimensions, float additionalParam, BSpline parametricCurve = BSpline(), bool in2D = false);
    static std::function<float(Vector3)> createPredefinedShapeFunction(PredefinedShapes shape, const Vector3& dimensions, float additionalParam, BSpline parametricCurve = BSpline(), bool in2D = false);
    static float isovalue;
    static float zResolution;


    static std::function<float(Vector3)> createSphereFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createBlockFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createGaussianFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createCylinderFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createRockFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createMountainFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createDuneFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createBasinFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createCaveFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createArchFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createNoise2DFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createMountainChainFunction(float sigma, float width, float depth, float height, BSpline path, bool in2D = false);
    static std::function<float(Vector3)> createPolygonFunction(float sigma, float width, float depth, float height, BSpline path, bool in2D = false);
    static std::function<float(Vector3)> createDistanceMapFunction(float sigma, float width, float depth, float height, BSpline path, bool in2D = false);
    static std::function<float(Vector3)> createParametricTunnelFunction(float sigma, float width, float depth, float height, BSpline path, bool in2D = false);
    static std::function<float(Vector3)> createRippleFunction(float sigma, float width, float depth, float height, bool in2D = false);
    static std::function<float(Vector3)> createIdentityFunction(float sigma, float width, float depth, float height, bool in2D = false);

//    static std::function<float(Vector3)> ...;

    static std::function<float(Vector3)> convert2DfunctionTo3Dfunction(std::function<float(Vector3)> func);

    ImplicitNaryOperator *getParent() const;
    void setParent(ImplicitNaryOperator* newParent);
    ImplicitNaryOperator* parent = nullptr;

    static int currentMaxIndex;

    static std::string json_identifier;

    virtual void reevaluateAll() {}

//protected:
    GridF _cachedMinHeight;
    GridF _cachedMaxHeight;
    GridF _cachedVoxelized;
    bool _cached = false;
};

class ImplicitPrimitive : public ImplicitPatch {
public:
    ImplicitPrimitive();

    float evaluate(const Vector3& pos) const;
    std::map<TerrainTypes, float> getMaterials(const Vector3& pos) const;

    AABBox getSupportBBox() const;
    AABBox getBBox() const;
    void update();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    void setDimensions(const Vector3& newDimensions);
    void setSupportDimensions(const Vector3& newSupportDimensions);

    bool contains(PredefinedShapes shape) { return this->predefinedShape == shape; }
    std::vector<ImplicitPatch*> findAll(PredefinedShapes shape);

    virtual ImplicitPatch* copy() const;

    Vector3 position = Vector3(false);
    Vector3 dimensions = Vector3(false);
    Vector3 supportDimensions = Vector3(false);
    std::function<float(Vector3)> evalFunction;
    std::function<float(Vector3)> evalFunction2D;
    TerrainTypes material = WATER;

    PredefinedShapes predefinedShape = None;
    std::vector<float> parametersProvided;

    std::string heightmapFilename = "";
    GridF cachedHeightmap;

    static ImplicitPrimitive* fromHeightmap(std::string filename, const Vector3& dimensions = Vector3(false), ImplicitPrimitive *prim = nullptr);
    static ImplicitPrimitive* fromHeightmap(GridF heightmap, std::string filename = "", ImplicitPrimitive *prim = nullptr);
};

class ImplicitNaryOperator : public ImplicitPatch {
public:
    ImplicitNaryOperator();

    float evaluate(const Vector3& pos) const;
    std::map<TerrainTypes, float> getMaterials(const Vector3& pos) const;
    AABBox getSupportBBox() const;
    AABBox getBBox() const;
    void update();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    bool contains(PredefinedShapes shape);
    virtual std::vector<ImplicitPatch*> findAll(PredefinedShapes shape);

    virtual void simplify();

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

    float evaluate(const Vector3& pos) const;
    std::map<TerrainTypes, float> getMaterials(const Vector3& pos) const;
    float evaluateFromAandB(float evalA, float evalB) const;

    float evaluateA(const Vector3& pos) const;
    float evaluateB(const Vector3& pos) const;

    std::map<TerrainTypes, float> getMaterialsA(const Vector3& pos) const;
    std::map<TerrainTypes, float> getMaterialsB(const Vector3& pos) const;
    std::pair<float, std::map<TerrainTypes, float>> getMaterialsAndTotalEvaluationA(const Vector3& pos) const;
    std::pair<float, std::map<TerrainTypes, float>> getMaterialsAndTotalEvaluationB(const Vector3& pos) const;

    bool contains(PredefinedShapes shape);
    virtual std::vector<ImplicitPatch*> findAll(PredefinedShapes shape);

    AABBox getSupportBBox() const;
    AABBox getBBox() const;
    void update();
    std::string toString();
    nlohmann::json toJson();
    static ImplicitPatch* fromJson(nlohmann::json content);

    void updateCache();

    void swapAB();


    Vector3 getEvaluationPositionForComposableA(const Vector3& pos) const;
    Vector3 getEvaluationPositionForComposableB(const Vector3& pos) const;

    virtual ImplicitPatch* copy() const;
    virtual void addChild(ImplicitPatch* newChild, int index);
    virtual void deleteAllChildren();

//    ImplicitPatch* composableA = nullptr;
//    ImplicitPatch* composableB = nullptr;
    ImplicitPatch* composableA() const;
    ImplicitPatch* composableB() const;

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

    float evaluate(const Vector3& point) const;

    void addTranslation(const Vector3& translation);

    void addRotation(const Vector3& rotationAngles);

    void addScaling(const Vector3& scaleFactors);

    void addDisplacementField(const GridV3& displacementField);

    Vector3 applyTransform(const Vector3& pos) const;
    Vector3 inverseTransform(const Vector3& pos) const;


    std::map<TerrainTypes, float> getMaterials(const Vector3& pos) const;

    AABBox getSupportBBox() const;
    AABBox getBBox() const;
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

    ImplicitPatch* composableA() const;

    void translate(const Vector3& translation);
    void rotate(const Vector3& angles); //float angleX, float angleY, float angleZ);
    void scale(const Vector3& scaleFactor);

    void addRandomNoise(float amplitude, float period = 20.f, float offset = 10.f);
    void addRandomWrap(float amplitude, float period = 20.f, float offset = 10.f);
    void addWrapFunction(GridV3 func);
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




class Implicit2DNary : public ImplicitNaryOperator {
public:
    Implicit2DNary();

    std::map<TerrainTypes, float> getMaterials(const Vector3& pos) const;
    std::pair<float, std::map<TerrainTypes, float>> getMaterialsAndTotalEvaluation(const Vector3 &pos) const;
    float evaluate(const Vector3& pos) const; // In this case, returns the height

    GridF getVoxelized(const Vector3& dimensions = Vector3(false), const Vector3& scale = Vector3(1.f, 1.f, 1.f));

    virtual bool contains(const Vector3& v);

    bool checkIsInGround(const Vector3& position);
    virtual Vector3 getNormal(const Vector3& pos) const;
//    virtual AABBox getSupportBBox() const;
//    virtual AABBox getBBox() const;
    virtual float getMaxHeight(const Vector3& pos);
    virtual float getMinHeight(const Vector3& pos);
    virtual float getMinimalHeight(const Vector3 &minBox = Vector3::min(), const Vector3& maxBox = Vector3::max());
    virtual float getMaximalHeight(const Vector3& minBox = Vector3::min(), const Vector3& maxBox = Vector3::max());

    void reevaluateAll();
    float computeHeight(const Vector3& pos) const;
};

class Implicit2DPrimitive : public ImplicitPrimitive {
public:
    Implicit2DPrimitive();

    float evaluate(const Vector3 &pos) const;
};





class UnaryOp {
public:
    UnaryOp();
    std::function<Vector3(Vector3)> wrap;
    std::function<Vector3(Vector3)> unwrap;
};
class UnaryOpTranslate: public UnaryOp {
public:
    UnaryOpTranslate(const Vector3& translation);
};
class UnaryOpRotate: public UnaryOp {
public:
    UnaryOpRotate(const Vector3& rotationAngles, const Vector3& center);
};
class UnaryOpScale: public UnaryOp {
public:
    UnaryOpScale(const Vector3& scaling, const Vector3& center);
};
class UnaryOpWrap: public UnaryOp {
public:
    UnaryOpWrap(FastNoiseLite noise, const Vector3& strength);
    UnaryOpWrap(std::function<Vector3(Vector3)> func);
    UnaryOpWrap(GridV3 wrapper);
    GridV3 wrapper;
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

    float evaluate(const Vector3& pos);
    std::map<TerrainMaterial, float> getMaterials(const Vector3& pos);
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
