#include "Octree.h"

OctreeNode::OctreeNode(const Vector3 &origin, const Vector3 &halfDimension)
    : origin(origin), halfDimension(halfDimension) {
    for (int i = 0; i < 8; ++i) children[i] = nullptr;
}

OctreeNode::~OctreeNode() {
    for (int i = 0; i < 8; ++i) delete children[i];
}

bool OctreeNode::intersects(const Vector3 &start, const Vector3 &end) const {
    // Check if the two AABBs are disjoint along any axis
    if (start.x > origin.x + halfDimension.x || end.x < origin.x - halfDimension.x) return false;
    if (start.y > origin.y + halfDimension.y || end.y < origin.y - halfDimension.y) return false;
    if (start.z > origin.z + halfDimension.z || end.z < origin.z - halfDimension.z) return false;

    // If the AABBs are not disjoint along any axis, they must intersect
    return true;
}

Octree::Octree() : root(nullptr)
{
}

Octree::Octree(const Vector3 &origin, const Vector3 &halfDimension)
    : root(new OctreeNode(origin, halfDimension))
{
//    std::cout << "Root : " << origin << " -- " << halfDimension << std::endl;
}

Octree::~Octree() {
    delete root;
}

void Octree::insert(const Vector3 &point, const int& pointIndex) {
    if (!Vector3::isInBox(point, this->root->origin - this->root->halfDimension, this->root->origin + this->root->halfDimension))
        return; // Don't add this point if it's not inside of the octree space
    insert(root, point, pointIndex);
}

std::vector<OctreeNodeData> Octree::queryRange(const Vector3 &start, const Vector3 &end) {
    Vector3 _start = Vector3::min(start, end);
    Vector3 _end = Vector3::max(start, end);
    std::vector<OctreeNodeData> result;
    queryRange(root, _start, _end, result);
    return result;
}

void Octree::insert(OctreeNode *node, const Vector3 &point, const int& pointIndex) {
    if (node->data.size() < 20 && node->children[0] == nullptr) {
        // If the node has no children and is not full, add the point here
        node->data.push_back(OctreeNodeData({point, pointIndex}));
    } else {
        // Otherwise, split the node and add the point to the appropriate child
        if (node->children[0] == nullptr) {
            for (int i = 0; i < 8; ++i) {
                Vector3 newOrigin = node->origin;
                newOrigin.x += node->halfDimension.x * ((i & 1) ? 0.5f : -0.5f);
                newOrigin.y += node->halfDimension.y * ((i & 2) ? 0.5f : -0.5f);
                newOrigin.z += node->halfDimension.z * ((i & 4) ? 0.5f : -0.5f);
                node->children[i] = new OctreeNode(newOrigin, node->halfDimension * 0.5f);
            }
        }

        auto dataCopy = node->data;
        node->data.clear();
        for (auto& data : dataCopy) {
            insert(node, data.pos, data.index);
        }

        // Add the point to the appropriate child
        int octant = (point.x >= node->origin.x) + ((point.y >= node->origin.y) << 1) + ((point.z >= node->origin.z) << 2);
        insert(node->children[octant], point, pointIndex);
    }
}

void Octree::queryRange(OctreeNode *node, const Vector3 &start, const Vector3 &end, std::vector<OctreeNodeData> &result) {
    // If the node does not intersect with the query range, return
    if (!node->intersects(start, end)) {
        return;
    }

    // If the node is a leaf node, check all points in the node
    if (node->children[0] == nullptr) {
        for (const auto& point : node->data) {
            if (Vector3::isInBox(point.pos, start, end)) {
                result.push_back(point);
            }
        }
    } else {
        // If the node is not a leaf node, check all children
        for (int i = 0; i < 8; ++i) {
            queryRange(node->children[i], start, end, result);
        }
    }
}
