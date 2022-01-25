#include "TreeColonisation.h"

using namespace TreeColonisationAlgo;

TreeColonisation::TreeColonisation()
{
}

TreeColonisation::TreeColonisation(std::vector<Vector3> nodes, Vector3 startPos, float segmentLength, float randomness)
{
    init(nodes, startPos, segmentLength, randomness);
}

TreeColonisation::TreeColonisation(Graph<NODE_TYPE> graph, Vector3 startPos, float segmentLength, float randomness)
{
    std::vector<Vector3> pos;
    for (const auto& node : graph.nodes)
        pos.push_back(node->pos);
    init(pos, startPos, segmentLength, randomness);
}

TreeColonisation::~TreeColonisation()
{
    this->segments.clear();
    this->nodes.clear();
}

void TreeColonisation::init(std::vector<Vector3> nodes, Vector3 startPos, float segmentLength, float randomness)
{
    this->segmentLength = segmentLength;
    this->treeSuccess = false;
    this->randomness = randomness;
    this->startPosition = startPos;
    this->initialNodes = nodes;
    Vector3 meanNodePos;
    this->nodes.clear();
    for (const auto& pos : nodes)
    {
        this->nodes.push_back(std::make_shared<Node>(pos, ATTRACTOR));
        meanNodePos += pos;
    }
    if (this->nodes.size() > 0)
        meanNodePos /= (float)this->nodes.size();

    // Make root with a random position but pointing towards the nodes...
    this->segments.clear();
    this->segments.push_back(std::make_shared<Segment>(nullptr, startPos, meanNodePos - startPos, this->segmentLength));
}

void TreeColonisation::reset(std::vector<Vector3> newKeypoints)
{
    if (!newKeypoints.empty()) {
        this->initialNodes.clear();
        this->initialNodes = newKeypoints;
    }
    this->init(this->initialNodes, this->startPosition, this->segmentLength, this->randomness);
}


void TreeColonisation::process()
{
    std::shared_ptr<Segment> current = this->segments[0];

    int maxIteration = 10000;
    this->treeSuccess = false;
    // Grow the trunk until we can see a nearby node
    while (current->distanceTo(this->nodes) > this->nodeMaxDistance)
    {
        current = current->next();
        this->segments.push_back(current);
        maxIteration --;
        if (maxIteration <= 0) {
            return;
        }
    }
    this->treeSuccess = true;

    maxIteration = 1000;
    size_t previousSegmentSize = 0;
    // Grow towards the leaves
    while (!this->nodes.empty())
    {
        // For each leaf, attract/repulse the closest segment.
        for (auto& node : this->nodes)
        {
            float closestDist = std::numeric_limits<float>::max();
            std::shared_ptr<Segment> closestSegment = nullptr;

            for (auto& segment : this->segments)
            {
                float returnedDistance = segment->distanceTo(node);

                // If the segment is "almost" touching the node, mark it as reached to be destroyed later
                // No need to continue
                if (returnedDistance < this->nodeMinDistance)
                {
                    node->reached = true;
                    break;
                }

                if (returnedDistance < closestDist) {
                    closestDist = returnedDistance;
                    closestSegment = segment;
                }
            }

            if (closestSegment != nullptr)
            {
                closestSegment->dir += (node->pos - closestSegment->pos).normalized() * (node->type == ATTRACTOR ? 1.f : -1.f);
                closestSegment->numberOfCurrentAttractors++;
            }
        }

        // Remove reached nodes
        for (int i = this->nodes.size() - 1; i >= 0; i--)
            if (this->nodes[i]->reached)
                this->nodes.erase(this->nodes.begin() + i);

        // Split branches if needed
        for (int i = this->segments.size() - 1; i >= 0; i--)
        {
            std::shared_ptr<Segment> seg = this->segments[i];
            if (seg->numberOfCurrentAttractors > 0)
            {
                seg->dir /= (float)seg->numberOfCurrentAttractors;
                seg->dir.normalize();
                if (seg->children.size() < 5)
                    this->segments.push_back(seg->next(this->randomness));
                seg->reset();
            }
        }
        maxIteration --;
        // If blocked in a loop (not able to find a new node to hang to), stop it
        if (previousSegmentSize == this->segments.size() || maxIteration <= 0) {
            this->nodes.clear();
            return;
        }
        previousSegmentSize = this->segments.size();
    }
    return;
}

#include <set>
std::vector<std::vector<Vector3> > TreeColonisation::simplifyPaths()
{
    std::vector<std::vector<Vector3> > allPaths;
    std::set<std::shared_ptr<Segment>> allSegments(this->segments.begin(), this->segments.end());
    std::set<std::shared_ptr<Segment>> Q;

    std::vector<Vector3> currentPath;
    if (this->treeSuccess)
    {
        Q.insert(this->segments[0]);
        while (!Q.empty())
        {
            std::shared_ptr<Segment> current = *Q.begin();
            Q.erase(Q.begin());
//            allSegments.erase(std::find(allSegments.begin(), allSegments.end(), current));
            if (current->parent != nullptr) {
                currentPath.push_back(current->parent->pos);
                currentPath.push_back(current->pos);
            }

            if (current->children.empty())
            {
                allPaths.push_back(currentPath);
                currentPath = std::vector<Vector3>();
            }
            for (auto& child : current->children)
            {
                Q.insert(Q.begin(), child);
//                allSegments.erase(std::find(allSegments.begin(), allSegments.end(), child));
            }

        }
    }
    return allPaths;
}

Segment::Segment()
{

}

Segment::Segment(std::shared_ptr<Segment> parent, Vector3 pos, Vector3 direction, float length)
    : parent(parent), pos(pos), dir(direction), originalDir(dir), length(length), numberOfCurrentAttractors(0)
{

}

float Segment::distanceTo(std::vector<std::shared_ptr<Node> > nodes)
{
    float minDist = std::numeric_limits<float>::max();
    for (auto& node : nodes) {
        float returnedDistance = this->distanceTo(node);
        if (returnedDistance < minDist) {
            minDist = returnedDistance;
        }
    }
    return minDist;
}

float Segment::distanceTo(std::shared_ptr<Node> node)
{
    return (node->pos - this->pos).norm();
}

std::shared_ptr<Segment> Segment::next(float randomness)
{
    Vector3 randomChildDir = Vector3(this->dir * (1-randomness) + Vector3::random()*randomness).normalized();

    std::shared_ptr<Segment> child = std::make_shared<Segment>(this->shared_from_this(), this->pos + this->dir * this->length, randomChildDir, this->length);
    this->children.push_back(child);
    return child;
}

void Segment::reset()
{
    this->dir = this->originalDir;
    this->numberOfCurrentAttractors = 0;
}

Node::Node()
{
}

Node::Node(Vector3 pos, NODE_TYPE type)
    : pos(pos), type(type), reached(false)
{

}
