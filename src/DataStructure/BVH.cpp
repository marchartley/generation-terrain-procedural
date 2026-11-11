#include "BVH.h"
#include "SpatialQueries.h"
#include <numeric>
#include <memory>
#include "Utils/Collisions.h"
#include "Utils/Utils.h"


BVHTree::BVHTree() : memoryPool(new BVHMemoryPool) {}

BVHTree::~BVHTree() { delete memoryPool; }

BVHTree::BVHTree(const BVHTree &other) : BVHTree() {
    maxTrianglesPerLeaves = other.maxTrianglesPerLeaves;
    useParallel = other.useParallel;
    useSAH = other.useSAH;
    useQuickSelect = other.useQuickSelect;
    try {
        build(other.triangles);
    } catch (...) {
        delete memoryPool;
        memoryPool = nullptr;
        throw;
    }
}

const BVHTree &BVHTree::operator=(const BVHTree &other) {
    if (this != &other) {
        BVHTree copy(other);
        swap(copy);
    }
    return *this;
}

void BVHTree::swap(BVHTree &other) noexcept {
    using std::swap;
    triangles.swap(other.triangles);
    order.swap(other.order);
    swap(root, other.root);
    swap(memoryPool, other.memoryPool);
    swap(maxTrianglesPerLeaves, other.maxTrianglesPerLeaves);
    swap(useParallel, other.useParallel);
    swap(useSAH, other.useSAH);
    swap(useQuickSelect, other.useQuickSelect);
}
BVHTree::BVHTree(BVHTree &&other) noexcept : memoryPool(nullptr) { swap(other); }
BVHTree &BVHTree::operator=(BVHTree &&other) noexcept {
    if (this != &other)
        swap(other);
    return *this;
}
void BVHTree::recordBuildFailure() noexcept {
#pragma omp critical(terrain_bvh_build_error)
    {
        if (!buildFailure)
            buildFailure = std::current_exception();
    }
}
const SpacePartitioning &BVHTree::build(const std::vector<Triangle> &input) {
    if (maxTrianglesPerLeaves < 1)
        throw std::invalid_argument("BVH leaf capacity must be positive");
    spatial_detail::validate(input);
    BVHTree replacement;
    replacement.maxTrianglesPerLeaves = maxTrianglesPerLeaves;
    replacement.useParallel = useParallel;
    replacement.useSAH = useSAH;
    replacement.useQuickSelect = useQuickSelect;
    replacement.triangles = input;
    replacement.order.resize(input.size());
    std::iota(replacement.order.begin(), replacement.order.end(), 0);
    if (!input.empty()) {
        if (useParallel && input.size() >= 1024) {
#pragma omp parallel
            {
#pragma omp single
                {
                    try {
                        replacement.root = replacement.buildBVH(0, int(input.size()));
                    } catch (...) {
                        replacement.recordBuildFailure();
                    }
                }
            }
        } else
            replacement.root = replacement.buildBVH(0, int(input.size()));
    }
    if (replacement.buildFailure)
        std::rethrow_exception(replacement.buildFailure);
    swap(replacement);
    return *this;
}

std::set<size_t> BVHTree::getAllStoredTrianglesIndices() const {
    return std::set<size_t>(order.begin(), order.end());
}

bool BVHTree::checkIntersection(const Vector3 &a, const Vector3 &b, const std::set<size_t> &ignored) const {
    spatial_detail::validate(a, b);
    return _checkIntersection(root, a, b, ignored);
}

bool BVHTree::_checkIntersection(BVHNode *node, const Vector3 &a, const Vector3 &b,
                                 const std::set<size_t> &ignored) const {
    double t;
    if (!node || !spatial_detail::box(node->box, a, b, 1, t))
        return false;
    for (auto i : node->trianglesIndices)
        if (!ignored.count(i) && spatial_detail::triangle(triangles[i], a, b, t))
            return true;
    return _checkIntersection(node->left, a, b, ignored) || _checkIntersection(node->right, a, b, ignored);
}

namespace {
AABBox rangeBox(const std::vector<Triangle> &ts, const std::vector<size_t> &order, int start, int end) {
    Vector3 lo = Vector3::max(), hi = Vector3::min();
    for (int i = start; i < end; ++i)
        for (const auto &v : ts[order[i]].vertices) {
            lo = Vector3::min(lo, v);
            hi = Vector3::max(hi, v);
        }
    return AABBox(lo, hi);
}
int largest(const AABBox &b) {
    auto d = b.dimensions();
    return d.x() > d.y() ? (d.x() > d.z() ? 0 : 2) : (d.y() > d.z() ? 1 : 2);
}
float centroid(const Triangle &t, int axis) {
    return float((double(t[0][axis]) + t[1][axis] + t[2][axis]) / 3);
}
double area(const AABBox &b) {
    auto d = b.dimensions();
    return 2 * (double(d.x()) * d.y() + double(d.x()) * d.z() + double(d.y()) * d.z());
}
} // namespace
int BVHTree::findBestSplitSAH(int start, int end) {
    if (start < 0 || end <= start || end - start < 2 || size_t(end) > order.size())
        throw std::out_of_range("BVH SAH range");
    int axis = largest(rangeBox(triangles, order, start, end)), n = end - start;
    std::sort(order.begin() + start, order.begin() + end, [&](size_t a, size_t b) {
        return centroid(triangles[a], axis) < centroid(triangles[b], axis);
    });
    std::vector<AABBox> right(n);
    right.back() = rangeBox(triangles, order, end - 1, end);
    for (int i = n - 2; i >= 0; --i) {
        right[i] = right[i + 1];
        for (auto &v : triangles[order[start + i]].vertices)
            right[i].expand(v);
    }
    AABBox left = rangeBox(triangles, order, start, start + 1);
    int split = start + 1;
    double cost = std::numeric_limits<double>::infinity();
    for (int i = 1; i < n; ++i) {
        double c = area(left) * i + area(right[i]) * (n - i);
        if (c < cost || (c == cost && std::abs(i - n / 2) < std::abs(split - start - n / 2))) {
            cost = c;
            split = start + i;
        }
        for (auto &v : triangles[order[start + i]].vertices)
            left.expand(v);
    }
    return split;
}

int BVHTree::partition(int start, int end, int pivot, int axis) {
    if (start < 0 || end <= start || size_t(end) > order.size() || pivot < start || pivot >= end ||
        axis < 0 || axis > 2)
        throw std::out_of_range("BVH partition range");
    auto value = centroid(triangles[order[pivot]], axis);
    return int(std::partition(order.begin() + start, order.begin() + end,
                              [&](size_t i) { return centroid(triangles[i], axis) < value; }) -
               order.begin());
}

// QuickSelect to find the median triangle along a given axis.
int BVHTree::quickSelect(int start, int end, int axis) {
    if (start < 0 || end < start || size_t(end) >= order.size() || axis < 0 || axis > 2)
        throw std::out_of_range("BVH partition range");
    int mid = start + (end - start + 1) / 2;
    std::nth_element(
        order.begin() + start, order.begin() + mid, order.begin() + end + 1,
        [&](size_t a, size_t b) { return centroid(triangles[a], axis) < centroid(triangles[b], axis); });
    return mid;
}

BVHNode *BVHTree::allocateNode() {
    if (!memoryPool)
        memoryPool = new BVHMemoryPool;
    return useParallel ? memoryPool->parallelAllocate() : memoryPool->allocate();
}

void BVHTree::traverseBVH(BVHNode *node, const Vector3 &p, float &distance, size_t &index,
                          const std::vector<Triangle> &ts) {
    if (!p.isValid() || std::isnan(distance) || distance < 0)
        throw std::invalid_argument("Invalid nearest query");
    double best = double(distance) * distance;
    auto visit = [&](auto &&self, BVHNode *current) -> void {
        if (!current || spatial_detail::distance2(current->box, p) > best)
            return;
        for (auto i : current->trianglesIndices) {
            double d = spatial_detail::triangleDistance2(ts.at(i), p);
            if (d < best || (d == best && i < index)) {
                best = d;
                index = i;
            }
        }
        auto a = current->left, b = current->right;
        if (a && b && spatial_detail::distance2(a->box, p) > spatial_detail::distance2(b->box, p))
            std::swap(a, b);
        self(self, a);
        self(self, b);
    };
    visit(visit, node);
    distance = float(std::sqrt(best));
}

size_t BVHTree::getClosestTriangle(const Vector3 &p) {
    if (!p.isValid())
        throw std::invalid_argument("Invalid nearest point");
    if (!root)
        throw std::out_of_range("Empty BVH");
    float distance = std::numeric_limits<float>::infinity();
    size_t index = size_t(-1);
    traverseBVH(root, p, distance, index, triangles);
    return index;
}


BVHNode *BVHTree::buildBVH(int start, int end) {
    if (start < 0 || end < start || size_t(end) > order.size() || maxTrianglesPerLeaves < 1)
        throw std::out_of_range("Invalid BVH build range/capacity");
    if (start == end)
        return nullptr;
    auto node = allocateNode();
    node->box = rangeBox(triangles, order, start, end);
    if (end - start <= maxTrianglesPerLeaves) {
        node->trianglesIndices.assign(order.begin() + start, order.begin() + end);
        return node;
    }
    int mid = useSAH ? findBestSplitSAH(start, end) : quickSelect(start, end - 1, largest(node->box));
    // Bound SAH recursion even for pathological distributions; median builds
    // already satisfy this balance bound.
    const int margin = std::max(1, (end - start) / 8);
    mid = std::clamp(mid, start + margin, end - margin);
    if (useParallel && end - start >= 1024) {
#pragma omp task shared(node) firstprivate(start, mid)
        {
            try {
                node->left = buildBVH(start, mid);
            } catch (...) {
                recordBuildFailure();
            }
        }
#pragma omp task shared(node) firstprivate(mid, end)
        {
            try {
                node->right = buildBVH(mid, end);
            } catch (...) {
                recordBuildFailure();
            }
        }
#pragma omp taskwait
    } else {
        node->left = buildBVH(start, mid);
        node->right = buildBVH(mid, end);
    }
    return node;
}

void BVHTree::closest(BVHNode *node, const Vector3 &a, const Vector3 &b, const std::set<size_t> &ignored,
                      double &best, std::pair<Vector3, size_t> &result) const {
    double entry;
    if (!node || !spatial_detail::box(node->box, a, b, best, entry))
        return;
    for (auto i : node->trianglesIndices) {
        if (ignored.count(i))
            continue;
        double t;
        if (spatial_detail::triangle(triangles[i], a, b, t) &&
            (t < best || (t == best && i < result.second))) {
            best = t;
            result = {spatial_detail::point(a, b, t), i};
        }
    }
    auto left = node->left, right = node->right;
    double l = 0, r = 0;
    bool hasLeft = left && spatial_detail::box(left->box, a, b, best, l),
         hasRight = right && spatial_detail::box(right->box, a, b, best, r);
    if (hasLeft && hasRight && r < l) {
        std::swap(left, right);
        std::swap(hasLeft, hasRight);
    }
    if (hasLeft)
        closest(left, a, b, ignored, best, result);
    if (hasRight)
        closest(right, a, b, ignored, best, result);
}
std::pair<Vector3, size_t> BVHTree::getIntersectionAndTriangleIndex(const Vector3 &a, const Vector3 &b,
                                                                    const std::set<size_t> &ignored) const {
    return _getIntersectionAndTriangleIndex(root, a, b, ignored);
}

std::pair<Vector3, size_t> BVHTree::_getIntersectionAndTriangleIndex(BVHNode *node, const Vector3 &a,
                                                                     const Vector3 &b,
                                                                     const std::set<size_t> &ignored) const {
    spatial_detail::validate(a, b);
    double best = 1;
    std::pair<Vector3, size_t> result{Vector3::invalid, size_t(-1)};
    closest(node, a, b, ignored, best, result);
    return result;
}

std::pair<Vector3, size_t> BVHTree::getIntersectionAndTriangleIndex(const Vector3 &a, const Vector3 &b,
                                                                    size_t ignored) const {
    return _getIntersectionAndTriangleIndex(root, a, b, ignored);
}

std::pair<Vector3, size_t> BVHTree::_getIntersectionAndTriangleIndex(BVHNode *node, const Vector3 &a,
                                                                     const Vector3 &b, size_t ignored) const {
    return _getIntersectionAndTriangleIndex(node, a, b, std::set<size_t>{ignored});
}

void BVHTree::all(BVHNode *node, const Vector3 &a, const Vector3 &b,
                  std::vector<std::pair<Vector3, size_t>> &result) const {
    double t;
    if (!node || !spatial_detail::box(node->box, a, b, 1, t))
        return;
    for (auto i : node->trianglesIndices)
        if (spatial_detail::triangle(triangles[i], a, b, t))
            result.emplace_back(spatial_detail::point(a, b, t), i);
    all(node->left, a, b, result);
    all(node->right, a, b, result);
}
std::vector<std::pair<Vector3, size_t>>
BVHTree::getAllIntersectionsAndTrianglesIndices(const Vector3 &a, const Vector3 &b) const {
    return _getAllIntersectionsAndTrianglesIndices(root, a, b);
}

std::vector<std::pair<Vector3, size_t>>
BVHTree::_getAllIntersectionsAndTrianglesIndices(BVHNode *node, const Vector3 &a, const Vector3 &b) const {
    spatial_detail::validate(a, b);
    std::vector<std::pair<Vector3, size_t>> result;
    all(node, a, b, result);
    return result;
}

BVHNode::BVHNode() : left(nullptr), right(nullptr) {}

BVHNode::~BVHNode() {
//    if (left != nullptr)
//        delete left;
//    if (right != nullptr)
//        delete right;
}
