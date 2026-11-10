#include "KDTree.h"
#include <algorithm>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>

KDNode::KDNode(size_t index, int dimension) : particleIndex(index), axis(dimension) {}
KDNode::~KDNode() {
    delete left;
    delete right;
}
namespace {
KDNode *clone(const KDNode *node) {
    if (!node)
        return nullptr;
    auto result = std::make_unique<KDNode>(node->particleIndex, node->axis);
    result->left = clone(node->left);
    result->right = clone(node->right);
    return result.release();
}
void checkQuery(const Vector3 &p, float radius) {
    if (!p.isValid() || !std::isfinite(radius) || radius < 0)
        throw std::invalid_argument("Invalid neighbor query");
}
} // namespace
KDTree::KDTree() = default;
KDTree::~KDTree() {
    delete root;
}
KDTree::KDTree(std::vector<Particle> &particles) {
    build(particles);
}
KDTree::KDTree(const KDTree &other) : positions(other.positions) {
    root = clone(other.root);
}
KDTree &KDTree::operator=(const KDTree &other) {
    if (this != &other) {
        KDTree copy(other);
        std::swap(root, copy.root);
        positions.swap(copy.positions);
    }
    return *this;
}
KDTree::KDTree(KDTree &&other) noexcept
    : root(std::exchange(other.root, nullptr)), positions(std::move(other.positions)) {}
KDTree &KDTree::operator=(KDTree &&other) noexcept {
    if (this != &other) {
        delete root;
        root = std::exchange(other.root, nullptr);
        positions = std::move(other.positions);
    }
    return *this;
}
KDNode *KDTree::build(std::vector<Particle> particles, int depth) {
    if (depth < 0)
        throw std::invalid_argument("Negative KDTree depth");
    std::vector<Vector3> coordinates;
    coordinates.reserve(particles.size());
    for (const auto &p : particles) {
        if (!p.position.isValid())
            throw std::invalid_argument("Invalid particle position");
        coordinates.push_back(p.position);
    }
    std::vector<size_t> indices(particles.size());
    std::iota(indices.begin(), indices.end(), 0);
    auto construct = [&](auto &&self, size_t begin, size_t end, int level) -> KDNode * {
        if (begin == end)
            return nullptr;
        size_t middle = begin + (end - begin) / 2;
        int axis = level % 3;
        std::nth_element(indices.begin() + begin, indices.begin() + middle, indices.begin() + end,
                         [&](size_t a, size_t b) { return coordinates[a][axis] < coordinates[b][axis]; });
        auto node = std::make_unique<KDNode>(indices[middle], axis);
        node->left = self(self, begin, middle, level + 1);
        node->right = self(self, middle + 1, end, level + 1);
        return node.release();
    };
    std::unique_ptr<KDNode> replacement(construct(construct, 0, indices.size(), depth));
    delete root;
    root = replacement.release();
    positions.swap(coordinates);
    return root;
}
void KDTree::query(KDNode *node, const Vector3 &p, double r2, std::vector<size_t> &result) const {
    if (!node)
        return;
    const auto &q = positions[node->particleIndex];
    double d = double(p[node->axis]) - q[node->axis];
    double dx = double(p.x()) - q.x(), dy = double(p.y()) - q.y(), dz = double(p.z()) - q.z();
    if (dx * dx + dy * dy + dz * dz <= r2)
        result.push_back(node->particleIndex);
    query(d < 0 ? node->left : node->right, p, r2, result);
    if (d * d <= r2)
        query(d < 0 ? node->right : node->left, p, r2, result);
}
std::vector<size_t> KDTree::findNeighbors(const Vector3 &p, float radius) const {
    checkQuery(p, radius);
    std::vector<size_t> result;
    query(root, p, double(radius) * radius, result);
    return result;
}
std::vector<size_t> KDTree::findNeighbors(std::vector<Particle> &particles, const Vector3 &p,
                                          float radius) const {
    if (particles.size() != positions.size())
        throw std::invalid_argument("Particle collection changed size; rebuild KDTree");
    return findNeighbors(p, radius);
}
void KDTree::findNeighbors(std::vector<Particle> &particles, KDNode *node, const Vector3 &p, float radius,
                           std::vector<size_t> &result) const {
    if (particles.size() != positions.size())
        throw std::invalid_argument("Particle collection changed size; rebuild KDTree");
    checkQuery(p, radius);
    query(node, p, double(radius) * radius, result);
}
