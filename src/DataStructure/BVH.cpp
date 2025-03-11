#include "BVH.h"
#include "Utils/Collisions.h"
#include "Utils/Utils.h"


BVHTree::BVHTree() : SpacePartitioning(), root(nullptr)
{
    memoryPool = new BVHMemoryPool();
}

BVHTree::~BVHTree()
{
    delete memoryPool;
//    if (root != nullptr)
//        delete root;
}

BVHTree::BVHTree(const BVHTree &other) : SpacePartitioning(other), root(nullptr) {
    memoryPool = new BVHMemoryPool(*other.memoryPool); // Deep copy
}

BVHTree &BVHTree::operator=(const BVHTree &other) {
    if (this != &other) { // Protect against self-assignment
        /*SpacePartitioning::operator=(other);*/ // call the base class's assignment operator

        // Clean up current resources
        /*delete memoryPool;*/

        // Deep copy
        /*root = other.root;*/
        // ... copy other members ...
        memoryPool = new BVHMemoryPool(/**other.memoryPool*/);
    }
    return *this;
}

SpacePartitioning& BVHTree::build(const std::vector<Triangle> &triangles)
{
    this->triangles = triangles;
    this->root = this->buildBVH(0, triangles.size());
    return *this;
}

std::set<size_t> BVHTree::getAllStoredTrianglesIndices() const
{
    std::vector<BVHNode*> queue = {root};
    std::vector<size_t> indices;

    while (!queue.empty()) {
        auto current = queue.back();
        queue.pop_back();
        indices.insert(indices.end(), current->trianglesIndices.begin(), current->trianglesIndices.end());
        if (current->left)
            queue.push_back(current->left);
        if (current->right)
            queue.push_back(current->right);
    }
    return convertVectorToSet(indices);
}

bool BVHTree::checkIntersection(const Vector3 &rayStart, const Vector3 &rayEnd, std::set<size_t> ignoredTriangles) const
{
    return this->_checkIntersection(root, rayStart, rayEnd, ignoredTriangles);
}

bool BVHTree::_checkIntersection(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd, std::set<size_t> ignoredTriangles) const
{
    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return false;
    }

    if (!node->trianglesIndices.empty()) { // leaf node
        for (const auto& triangleIndex : node->trianglesIndices) {
            if (isIn(triangleIndex, ignoredTriangles)) continue;
            auto& triangle = this->triangles[triangleIndex];
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid()) {
                return true;
            }
        }
    } else { // non-leaf node
        if (_checkIntersection(node->left, rayStart, rayEnd, ignoredTriangles))
            return true;
        if (_checkIntersection(node->right, rayStart, rayEnd, ignoredTriangles))
            return true;
    }

    return false;
}

float computeSurfaceArea(const AABBox& box) {
    Vector3 dims = box.dimensions();
    return 2.0f * (dims.x * dims.y + dims.x * dims.z + dims.y * dims.z);
}
int BVHTree::findBestSplitSAH(int start, int end) {
    const int MAX_SPLITS = 100; // maximum number of splits to consider
    int triangleCount = end - start;
    int step = std::max(1, triangleCount / MAX_SPLITS);

    std::vector<AABBox> leftCumulativeBoxes(end - start);
    std::vector<AABBox> rightCumulativeBoxes(end - start);

    // Compute left cumulative bounding boxes
    leftCumulativeBoxes[0].expand(triangles[start].vertices);
    for (int i = start + 1; i < end; ++i) {
        leftCumulativeBoxes[i - start] = leftCumulativeBoxes[i - start - 1];
        leftCumulativeBoxes[i - start].expand(triangles[i].vertices);
    }

    // Compute right cumulative bounding boxes
    rightCumulativeBoxes[end - start - 1].expand(triangles[end - 1].vertices);
    for (int i = end - 2; i >= start; --i) {
        rightCumulativeBoxes[i - start] = rightCumulativeBoxes[i - start + 1];
        rightCumulativeBoxes[i - start].expand(triangles[i].vertices);
    }

    float bestCost = std::numeric_limits<float>::max();
    int bestSplit = start;

    // Declare shared variables to hold the best cost and best split across threads
    float globalBestCost = bestCost;
    int globalBestSplit = bestSplit;

    #pragma omp parallel for shared(globalBestCost, globalBestSplit)
    for (int i = start; i < end; i += step) {
        AABBox leftBox = leftCumulativeBoxes[i - start];
        AABBox rightBox = (i + 1 < end) ? rightCumulativeBoxes[i - start + 1] : AABBox();  // If i+1 is out of bounds, use an empty box

        float leftSA = computeSurfaceArea(leftBox);
        float rightSA = computeSurfaceArea(rightBox);
        float cost = leftSA * (i - start + 1) + rightSA * (end - i - 1);

        // Critical section to update the global best cost and best split if the current thread has a better value
        #pragma omp critical
        {
            if (cost < globalBestCost) {
                globalBestCost = cost;
                globalBestSplit = i;
            }
        }
    }

    return globalBestSplit;
}

int BVHTree::partition(int start, int end, int pivotIdx, int axis) {
    Triangle& pivot = triangles[pivotIdx];
    Vector3 pivotMidPoint = (pivot[0] + pivot[1] + pivot[2]) / 3;

    // Move the pivot value to the end
    std::swap(triangles[pivotIdx], triangles[end-1]);

    int storeIndex = start;

    for (int i = start; i < end-1; ++i) {
        Vector3 triangleMidPoint = (triangles[i][0] + triangles[i][1] + triangles[i][2]) / 3;
        if (triangleMidPoint[axis] < pivotMidPoint[axis]) {
            std::swap(triangles[i], triangles[storeIndex]);
            storeIndex++;
        }
    }

    // Move pivot to its final place
    std::swap(triangles[storeIndex], triangles[end-1]);

    return storeIndex;
}

// QuickSelect to find the median triangle along a given axis.
int BVHTree::quickSelect(int start, int end, int axis) {
    if (start == end) {
        return start;
    }

    // Choose pivot randomly
    int pivotIdx = start + rand() % (end - start + 1);
    pivotIdx = partition(start, end, pivotIdx, axis);

    int median = start + (end - start) / 2;

    if (median == pivotIdx) {
        return median;
    } else if (median < pivotIdx) {
        return quickSelect(start, pivotIdx - 1, axis);
    } else {
        return quickSelect(pivotIdx + 1, end, axis);
    }
}

BVHNode *BVHTree::allocateNode() {
    if (useParallel) {
        return this->memoryPool->parallelAllocate();
    } else {
        return this->memoryPool->allocate();
    }
}

void BVHTree::traverseBVH(BVHNode *node, const Vector3 &pos, float &minDistance, size_t &closestTriangleIndex, const std::vector<Triangle> &triangles) {
    // If the current node is a leaf node
    if (!node->left && !node->right) {
        // Check distance to each triangle in the leaf node
        for (size_t triangleIndex : node->trianglesIndices) {
//            float dist = distanceToTriangle(pos, triangles[triangleIndex]);
            const Triangle tri = triangles[triangleIndex];
            float dist = std::sqrt(Collision::pointToTriangleDistanceSquared(tri[0], tri[1], tri[2], pos));
            if (dist < minDistance) {
                minDistance = dist;
                closestTriangleIndex = triangleIndex;
            }
        }
    } else {
        // Check if the point is closer to the left or right child node
        float leftDist = node->left->box.distanceTo(pos);
        float rightDist = node->right->box.distanceTo(pos);
        BVHNode* closerChild = (leftDist < rightDist) ? node->left : node->right;
        BVHNode* fartherChild = (leftDist < rightDist) ? node->right : node->left;

        // Recursively traverse the closer child node first
        traverseBVH(closerChild, pos, minDistance, closestTriangleIndex, triangles);

        // Check if the farther child node needs to be visited
        if (rightDist < minDistance) {
            traverseBVH(fartherChild, pos, minDistance, closestTriangleIndex, triangles);
        }
    }
}

size_t BVHTree::getClosestTriangle(const Vector3 &pos/*, const std::vector<Triangle> &triangles*/) {
    // Initialize minimum distance and closest triangle index
    float minDistance = std::numeric_limits<float>::max();
    size_t closestTriangleIndex = 0;

    // Start traversal from the root node
    traverseBVH(root, pos, minDistance, closestTriangleIndex, triangles);

    return closestTriangleIndex;
}


BVHNode *BVHTree::buildBVH(int start, int end)
{
    BVHNode* node = allocateNode();

    // Compute bounding box of all triangles from start to end
    Vector3 minPoint = Vector3::max();
    Vector3 maxPoint = Vector3::min();
    for (int i = start; i < end; ++i) {
        for (int ii = 0; ii < 3; ii++) {
            minPoint = Vector3::min(minPoint, triangles[i][ii]);
            maxPoint = Vector3::max(maxPoint, triangles[i][ii]);
        }
    }
    node->box = AABBox(minPoint, maxPoint);

    // Base case: if number of triangles is less than or equal to maxTrianglesPerLeaves, make a leaf node
    if (end - start <= maxTrianglesPerLeaves) {
        node->trianglesIndices = std::vector<size_t>(end - start);
        for (int i = 0; i < (end - start); i++)
            node->trianglesIndices[i] = start + i;
    } else {
        // Choose an axis and midpoint along that axis to partition the triangles
        Vector3 boxSize = maxPoint - minPoint;
        int axis = boxSize.x > boxSize.y ? (boxSize.x > boxSize.z ? 0 : 2) : (boxSize.y > boxSize.z ? 1 : 2);

        int mid;
        if (useSAH) {
            // Sort the triangles based on their midpoint along the chosen axis
            std::sort(triangles.begin() + start, triangles.begin() + end,
                [axis](const Triangle& triangle1, const Triangle& triangle2) {
                    Vector3 midPoint1 = (triangle1[0] + triangle1[1] + triangle1[2]) / 3;
                    Vector3 midPoint2 = (triangle2[0] + triangle2[1] + triangle2[2]) / 3;
                    return midPoint1[axis] < midPoint2[axis];
                });
            mid = this->findBestSplitSAH(start, end);
            if (mid == start || mid == end) {
                // Create a leaf node if we can't split further
                node->trianglesIndices = std::vector<size_t>(end - start);
                for (int i = 0; i < (end - start); i++)
                    node->trianglesIndices[i] = start + i;
            } else {
                node->left = buildBVH(start, mid);
                node->right = buildBVH(mid, end);
            }
            return node;
        } else if (useQuickSelect) {
            mid = this->quickSelect(start, end-1, axis);
            if (mid == start || mid == end - 1) {
                mid = start + (end - start) / 2;
            }
        } else {
            // Sort the triangles based on their midpoint along the chosen axis
            std::sort(triangles.begin() + start, triangles.begin() + end,
                [axis](const Triangle& triangle1, const Triangle& triangle2) {
                    Vector3 midPoint1 = (triangle1[0] + triangle1[1] + triangle1[2]) / 3;
                    Vector3 midPoint2 = (triangle2[0] + triangle2[1] + triangle2[2]) / 3;
                    return midPoint1[axis] < midPoint2[axis];
                });
            mid = start + (end - start) / 2;
        }

        if (useParallel) {
            // Parallelize the BVH node construction for the left and right children using OpenMP
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    node->left = buildBVH(start, mid);
                }
                #pragma omp section
                {
                    node->right = buildBVH(mid, end);
                }
            }
        } else {
            node->left = buildBVH(start, mid);
            node->right = buildBVH(mid, end);
        }
    }

    return node;
}

std::pair<Vector3, size_t> BVHTree::getIntersectionAndTriangleIndex(const Vector3 &rayStart, const Vector3 &rayEnd, std::set<size_t> ignoredTriangles) const
{
    return this->_getIntersectionAndTriangleIndex(root, rayStart, rayEnd, ignoredTriangles);
}

std::pair<Vector3, size_t> BVHTree::_getIntersectionAndTriangleIndex(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd, std::set<size_t> ignoredTriangles) const
{
    std::pair<Vector3, size_t> result = {Vector3(false), -1};
    float closestDistance = std::numeric_limits<float>::max();

    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return result; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->trianglesIndices.empty()) { // leaf node
        for (const auto& triangleIndex : node->trianglesIndices) {
            if (isIn(triangleIndex, ignoredTriangles)) continue;
            auto& triangle = this->triangles[triangleIndex];
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], false);
            if (intersectionPoint.isValid() && (intersectionPoint - rayStart).norm2() < closestDistance) {
                result = {intersectionPoint, triangleIndex};
                closestDistance = (intersectionPoint - rayStart).norm2();
            }
        }
    } else { // non-leaf node
        auto leftIntersection = _getIntersectionAndTriangleIndex(node->left, rayStart, rayEnd, ignoredTriangles);
        auto rightIntersection = _getIntersectionAndTriangleIndex(node->right, rayStart, rayEnd, ignoredTriangles);
        if (!leftIntersection.first.isValid()) {
            result = rightIntersection;
        } else if (!rightIntersection.first.isValid()) {
            result = leftIntersection;
        } else if ((leftIntersection.first - rayStart) < (rightIntersection.first - rayStart)) {
            result = leftIntersection;
        } else {
            result = rightIntersection;
        }
    }

    return result;
}

std::pair<Vector3, size_t> BVHTree::getIntersectionAndTriangleIndex(const Vector3 &rayStart, const Vector3 &rayEnd, size_t ignoredTriangle) const
{
    return this->_getIntersectionAndTriangleIndex(this->root, rayStart, rayEnd, ignoredTriangle);
}

std::pair<Vector3, size_t> BVHTree::_getIntersectionAndTriangleIndex(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd, size_t ignoredTriangle) const
{
    std::pair<Vector3, size_t> result = {Vector3(false), -1};
    float closestDistance = std::numeric_limits<float>::max();

    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return result; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->trianglesIndices.empty()) { // leaf node
        for (const auto& triangleIndex : node->trianglesIndices) {
            if (triangleIndex == ignoredTriangle) continue;
            auto& triangle = this->triangles[triangleIndex];
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid() && (intersectionPoint - rayStart).norm2() < closestDistance) {
                result = {intersectionPoint, triangleIndex};
                closestDistance = (intersectionPoint - rayStart).norm2();
            }
        }
    } else { // non-leaf node
        auto leftIntersection = _getIntersectionAndTriangleIndex(node->left, rayStart, rayEnd, ignoredTriangle);
        auto rightIntersection = _getIntersectionAndTriangleIndex(node->right, rayStart, rayEnd, ignoredTriangle);
        if (!leftIntersection.first.isValid()) {
            result = rightIntersection;
        } else if (!rightIntersection.first.isValid()) {
            result = leftIntersection;
        } else if ((leftIntersection.first - rayStart) < (rightIntersection.first - rayStart)) {
            result = leftIntersection;
        } else {
            result = rightIntersection;
        }
    }

    return result;
}

std::vector<std::pair<Vector3, size_t> > BVHTree::getAllIntersectionsAndTrianglesIndices(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    return this->_getAllIntersectionsAndTrianglesIndices(root, rayStart, rayEnd);
}

std::vector<std::pair<Vector3, size_t> > BVHTree::_getAllIntersectionsAndTrianglesIndices(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    std::vector<std::pair<Vector3, size_t>> result;
    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return result; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->trianglesIndices.empty()) { // leaf node
        for (const auto& triangleIndex : node->trianglesIndices) {
            auto& triangle = this->triangles[triangleIndex];
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid()) {
                result.push_back({intersectionPoint, triangleIndex});
//                closestDistance = (intersectionPoint - rayStart).norm2();
            }
        }
    } else { // non-leaf node
        auto leftIntersection = _getAllIntersectionsAndTrianglesIndices(node->left, rayStart, rayEnd);
        auto rightIntersection = _getAllIntersectionsAndTrianglesIndices(node->right, rayStart, rayEnd);

        result.insert(result.end(), leftIntersection.begin(), leftIntersection.end());
        result.insert(result.end(), rightIntersection.begin(), rightIntersection.end());
    }

    return result;
}

BVHNode::BVHNode() : left(nullptr), right(nullptr) {}

BVHNode::~BVHNode() {
//    if (left != nullptr)
//        delete left;
//    if (right != nullptr)
//        delete right;
}
/*
BVHMemoryPool::BVHMemoryPool() {
    allocateBlock();
}

BVHMemoryPool::~BVHMemoryPool() {
    for (NodeBlock* block : blocks) {
        delete block;  // Delete each allocated block
    }
}

BVHNode* BVHMemoryPool::allocateNode() {
    if (currentBlock->usedNodes == NodeBlock::BLOCK_SIZE) {
        allocateBlock();
    }
    return &(currentBlock->nodes[currentBlock->usedNodes++]);
}

BVHNode *BVHMemoryPool::parallelAllocateNode()
{
    BVHNode* res;
    #pragma omp critical
    {
        if (currentBlock->usedNodes == NodeBlock::BLOCK_SIZE) {
            allocateBlock();
        }
        res = &(currentBlock->nodes[currentBlock->usedNodes++]);
    }
    return res;
}

void BVHMemoryPool::allocateBlock() {
    currentBlock = new NodeBlock();  // Use new to allocate
    blocks.emplace_back(currentBlock);
}

void BVHMemoryPool::parallelAllocateBlock()
{
    #pragma omp critical
    {
        currentBlock = new NodeBlock();  // Use new to allocate
        blocks.emplace_back(currentBlock);
    }
}

NodeBlock::NodeBlock()
{
    this->nodes = new BVHNode[NodeBlock::BLOCK_SIZE];
}

NodeBlock::~NodeBlock()
{
    delete[] nodes;
}
*/
