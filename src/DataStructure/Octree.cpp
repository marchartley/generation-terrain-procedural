#include "Octree.h"
#include "SpatialQueries.h"
#include <memory>
#include <numeric>
#include <utility>
#include "Utils/Collisions.h"
#include "Utils/Utils.h"
#include <set>

OctreeNode::OctreeNode(const Vector3 &center, const Vector3 &half)
    : origin(center), halfDimension(half), sphereSqrRadius(half.norm2()) {
    std::fill(std::begin(children), std::end(children), nullptr);
}

OctreeNode::~OctreeNode() {
    for (auto child : children)
        delete child;
}

bool OctreeNode::intersects(const Vector3 &start, const Vector3 &end) const {
    auto lo = origin - halfDimension, hi = origin + halfDimension;
    for (int i = 0; i < 3; ++i)
        if (start[i] > hi[i] || end[i] < lo[i])
            return false;
    return true;
}

bool OctreeNode::intersects(const Vector3 &a, const Vector3 &b, const Vector3 &c) const {
    return intersects(Vector3::min({a, b, c}), Vector3::max({a, b, c}));
}

bool OctreeNode::intersects(const Triangle &t) const {
    return intersects(t[0], t[1], t[2]);
}

bool OctreeNode::sphereIntersectsTriangle(const Triangle &t) const {
    return spatial_detail::triangleDistance2(t, origin) <= sphereSqrRadius;
}

namespace {
std::unique_ptr<OctreeNode> clone(const OctreeNode *node) {
    if (!node)
        return {};
    auto result = std::make_unique<OctreeNode>(node->origin, node->halfDimension);
    result->data = node->data;
    for (int i = 0; i < 8; ++i)
        result->children[i] = clone(node->children[i]).release();
    return result;
}
Vector3 childCenter(const OctreeNode &n, int i) {
    return n.origin +
           n.halfDimension * Vector3((i & 1) ? .5f : -.5f, (i & 2) ? .5f : -.5f, (i & 4) ? .5f : -.5f);
}
std::unique_ptr<OctreeNode> construct(const Vector3 &center, const Vector3 &half,
                                      const std::vector<Triangle> &ts, const std::vector<size_t> &ids,
                                      int capacity, int depth) {
    auto node = std::make_unique<OctreeNode>(center, half);
    if (ids.size() <= size_t(capacity) || depth >= 24 ||
        half.maxComp() <= std::numeric_limits<float>::epsilon() * std::max(1.f, center.abs().maxComp())) {
        for (auto i : ids)
            node->data.emplace_back(ts[i], int(i));
        return node;
    }
    std::vector<size_t> bins[8];
    for (auto i : ids) {
        const auto &t = ts[i];
        int selected = -1;
        for (int c = 0; c < 8; ++c) {
            Vector3 mid = childCenter(*node, c), h = half * .5f;
            AABBox bounds(mid - h, mid + h);
            if (bounds.contains(t[0]) && bounds.contains(t[1]) && bounds.contains(t[2])) {
                selected = c;
                break;
            }
        }
        // Straddlers stay in their parent: no duplication or explosive subdivision.
        if (selected < 0)
            node->data.emplace_back(t, int(i));
        else
            bins[selected].push_back(i);
    }
    for (int c = 0; c < 8; ++c)
        if (!bins[c].empty())
            node->children[c] =
                construct(childCenter(*node, c), half * .5f, ts, bins[c], capacity, depth + 1).release();
    return node;
}
void candidates(const OctreeNode *node, const Vector3 &a, const Vector3 &b, std::vector<size_t> &result) {
    if (!node)
        return;
    double entry;
    if (!spatial_detail::box(AABBox(node->origin - node->halfDimension, node->origin + node->halfDimension),
                             a, b, 1, entry))
        return;
    for (const auto &d : node->data)
        result.push_back(size_t(d.index));
    for (auto child : node->children)
        candidates(child, a, b, result);
}
} // namespace
Octree::Octree() : SpacePartitioning(), root(nullptr)
{
}

Octree::Octree(const Vector3 &origin, const Vector3 &half, int capacity) : maxDataCapacity(capacity) {
    if (!origin.isValid() || !half.isValid() || half.minComp() < 0 || capacity < 1)
        throw std::invalid_argument("Invalid octree bounds/capacity");
    root = new OctreeNode(origin, half);
}


Octree::Octree(const Octree &other) : SpacePartitioning(other), maxDataCapacity(other.maxDataCapacity) {
    root = clone(other.root).release();
}
Octree &Octree::operator=(const Octree &other) {
    if (this != &other) {
        Octree copy(other);
        triangles.swap(copy.triangles);
        std::swap(root, copy.root);
        std::swap(maxDataCapacity, copy.maxDataCapacity);
    }
    return *this;
}
Octree::Octree(Octree &&other) noexcept
    : root(std::exchange(other.root, nullptr)), maxDataCapacity(other.maxDataCapacity) {
    triangles = std::move(other.triangles);
}
Octree &Octree::operator=(Octree &&other) noexcept {
    if (this != &other) {
        delete root;
        root = std::exchange(other.root, nullptr);
        triangles = std::move(other.triangles);
        maxDataCapacity = other.maxDataCapacity;
    }
    return *this;
}
Octree::~Octree() {
    delete root;
}

bool Octree::insert(std::vector<Triangle> added) {
    if (added.empty())
        return false;
    auto all = triangles;
    all.insert(all.end(), added.begin(), added.end());
    build(all);
    return true;
}

SpacePartitioning &Octree::build(const std::vector<Triangle> &input) {
    if (maxDataCapacity < 1)
        throw std::invalid_argument("Octree capacity must be positive");
    spatial_detail::validate(input);
    std::vector<Triangle> copy(input);
    std::unique_ptr<OctreeNode> replacement;
    if (!input.empty()) {
        Vector3 lo = Vector3::max(), hi = Vector3::min();
        for (auto &t : input)
            for (auto &v : t.vertices) {
                lo = Vector3::min(lo, v);
                hi = Vector3::max(hi, v);
            }
        std::vector<size_t> ids(input.size());
        std::iota(ids.begin(), ids.end(), 0);
        Vector3 center((double(lo.x()) + hi.x()) * .5, (double(lo.y()) + hi.y()) * .5,
                       (double(lo.z()) + hi.z()) * .5);
        Vector3 half = Vector3::max(center - lo, hi - center);
        replacement = construct(center, half, input, ids, maxDataCapacity, 0);
    }
    delete root;
    root = replacement.release();
    triangles.swap(copy);
    return *this;
}

std::set<size_t> Octree::getAllStoredTrianglesIndices() const {
    std::set<size_t> result;
    std::vector<OctreeNode *> stack;
    if (root)
        stack.push_back(root);
    while (!stack.empty()) {
        auto node = stack.back();
        stack.pop_back();
        for (auto &d : node->data)
            result.insert(d.index);
        for (auto child : node->children)
            if (child)
                stack.push_back(child);
    }
    return result;
}

std::pair<Vector3, size_t> Octree::getIntersectionAndTriangleIndex(const Vector3 &a, const Vector3 &b,
                                                                   const std::set<size_t> &ignored) const {
    spatial_detail::validate(a, b);
    std::vector<size_t> ids;
    candidates(root, a, b, ids);
    double best = std::numeric_limits<double>::infinity();
    std::pair<Vector3, size_t> result{Vector3::invalid, size_t(-1)};
    for (auto i : ids) {
        if (ignored.count(i))
            continue;
        double t;
        if (spatial_detail::triangle(triangles.at(i), a, b, t) &&
            (t < best || (t == best && i < result.second))) {
            best = t;
            result = {spatial_detail::point(a, b, t), i};
        }
    }
    return result;
}

std::pair<Vector3, size_t> Octree::getIntersectionAndTriangleIndex(const Vector3 &a, const Vector3 &b,
                                                                   size_t ignored) const {
    return getIntersectionAndTriangleIndex(a, b, std::set<size_t>{ignored});
}

std::vector<std::pair<Vector3, size_t>>
Octree::getAllIntersectionsAndTrianglesIndices(const Vector3 &a, const Vector3 &b) const {
    spatial_detail::validate(a, b);
    std::vector<size_t> ids;
    candidates(root, a, b, ids);
    std::vector<std::pair<Vector3, size_t>> result;
    for (auto i : ids) {
        double t;
        if (spatial_detail::triangle(triangles.at(i), a, b, t))
            result.emplace_back(spatial_detail::point(a, b, t), i);
    }
    return result;
}

std::vector<OctreeNodeData> Octree::queryRange(const Vector3 &start, const Vector3 &end) const {
    spatial_detail::validate(start, end);
    std::vector<OctreeNodeData> result;
    queryRange(root, Vector3::min(start, end), Vector3::max(start, end), result);
    return result;
}



bool Octree::insert(OctreeNode *node, const Vector3 &a, const Vector3 &b, const Vector3 &c,
                    const int &index) {
    if (!node || !node->intersects(a, b, c))
        return false;
    if (index < 0 || size_t(index) >= triangles.size())
        throw std::out_of_range("Octree triangle index");
    node->data.emplace_back(a, b, c, index);
    return true;
}

void Octree::queryRange(OctreeNode *node, const Vector3 &start, const Vector3 &end,
                        std::vector<OctreeNodeData> &result, int) const {
    if (!node || !node->intersects(start, end))
        return;
    result.insert(result.end(), node->data.begin(), node->data.end());
    for (auto child : node->children)
        queryRange(child, start, end, result);
}
