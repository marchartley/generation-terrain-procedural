#ifndef TREECOLONISATION_H
#define TREECOLONISATION_H

#include <vector>
#include "DataStructure/Vector3.h"
#include "Graph/FastPoissonGraph.h"

namespace TreeColonisationAlgo {

    enum NODE_TYPE
    {
        ATTRACTOR = 0,
        REPULSOR = 1
    };

    class Node;
    class TreeColonisation;
    class Segment;

    class TreeColonisation
    {
    public:
        TreeColonisation();
        TreeColonisation(std::vector<Vector3> nodes, const Vector3& startPos, float segmentLength = 1.f, float randomness = 0.f);
        TreeColonisation(GraphTemplate<NODE_TYPE> graph, const Vector3& startPos, float segmentLength = 1.f, float randomness = 0.f);
        ~TreeColonisation();

        void init(std::vector<Vector3> nodes, const Vector3& startPos, float segmentLength, float randomness);
        void reset(std::vector<Vector3> newKeypoints = std::vector<Vector3>());

        void process();

        std::pair<std::string, std::string> toFile();

        std::vector<std::vector<Vector3>> simplifyPaths();

        float nodeMinDistance;
        float nodeMaxDistance;

        Vector3 startPosition;
        float segmentLength;
        float randomness;
        std::vector<Vector3> initialNodes;

        std::vector<std::shared_ptr<Node>> nodes;
        std::vector<std::shared_ptr<Segment>> segments;

        bool treeSuccess;
    };

    class Segment : public std::enable_shared_from_this<Segment> {
    public:
        Segment();
        Segment(std::shared_ptr<Segment> parent, const Vector3& pos, const Vector3& direction, float length);

        float distanceTo(std::vector<std::shared_ptr<Node>> nodes);
        float distanceTo(std::shared_ptr<Node> node);

        std::shared_ptr<Segment> next(float randomness = 0.f);

        void reset();

        std::shared_ptr<Segment> parent;
        Vector3 pos;
        Vector3 dir;
        Vector3 originalDir;
        float length;

        int numberOfCurrentAttractors;
        std::vector<std::shared_ptr<Segment>> children;
    };

    class Node {
    public:
        Node();
        Node(const Vector3& pos, NODE_TYPE type = ATTRACTOR);

        Vector3 pos;
        NODE_TYPE type;
        bool reached;
    };

}
#endif // TREECOLONISATION_H
