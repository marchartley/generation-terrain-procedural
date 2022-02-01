#ifndef KARSTPATHSGENERATION_H
#define KARSTPATHSGENERATION_H

#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"
#include "Graph/FastPoissonGraph.h"

enum KARST_NODE_TYPE
{
    NORMAL = 0,
    FORCED_POINT = 1,
    DEAD_END = 2,
    SOURCE = 3,
    EXIT = 4
};

struct WaterHeight {
    WaterHeight(float height, float max_dist = 50.f, float effect = 1.f);
    float cost(Vector3 node) const;
    float height;
    float effect = 1.f;
    float max_distance_of_effect = 50.f;
};
struct PorositySphere {
    PorositySphere(Vector3 pos, float radius, float effect = 1.f);
    float cost(Vector3 node) const;
    Vector3 pos;
    float radius;
    float effect = 1.f;
};
struct FractureDirection {
    FractureDirection(Vector3 direction, float effect = 1.f);
    float cost(Vector3 nodeA, Vector3 nodeB) const;
    Vector3 direction;
    float effect = 1.f;
};

class KarstPathsGeneration
{
public:
    KarstPathsGeneration();
    KarstPathsGeneration(Matrix3<int> availablityMap, Vector3 karst_dimensions, float poissonDistance = 1.f, float gamma = 1.f);

    void resetSpecialNodes();
    void setSpecialNode(int node, KARST_NODE_TYPE mode);

    float computeCost(int nodeA, int nodeB);

    void addWaterHeight(WaterHeight waterHeight);
    void addPorositySphere(PorositySphere porositySphere);
    void addFractureDirection(FractureDirection fractureDirection);
    void addPorosityMap(Matrix3<float> porosity);

    void createEdges(int maxNumberNeighbors, float maxNeighboringDistance);

    void computeAllPathsBetweenSpecialNodes(int uniqueNodeToRecompute = -1);

    void updateTortuosity(float distanceToDisplace, std::vector<int> ignoreForSpecificNodes = std::vector<int>());

    Vector3 getNodePos(int node);

    static float smoothStep(float x) { return 3*std::pow(x, 2) - 2*std::pow(x, 3); }

    FastPoissonGraph<int> graph;

    std::vector<Vector3> tortuosityOffsets;

    std::vector<std::tuple<int, KARST_NODE_TYPE>> specialNodes;
    std::vector<int> mainPath;
    Matrix3<std::vector<int>> allPaths;
    std::vector<std::vector<int>> finalPaths;

    Matrix3<int> karst_available_matrix;

    float gamma;

    float distanceWeight, fractureWeight, porosityWeight, waterHeightWeight;
    std::vector<WaterHeight> waterHeights;
    std::vector<PorositySphere> porositySpheres;
    std::vector<FractureDirection> fracturesDirections;
    Matrix3<float> porosityMap;
};

#endif // KARSTPATHSGENERATION_H
