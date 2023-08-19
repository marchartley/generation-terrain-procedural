#include "KarstHole.h"
#include "Utils/Collisions.h"
#include "Utils/Utils.h"

KarstHole::KarstHole(float width, float height, KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
{
    this->startingProfile = KarstHoleProfile(startingShape).setSize(width, height);
    this->endingProfile = KarstHoleProfile(endingShape).setSize(width, height);
    this->path = BSpline();
    this->width = width;
    this->height = height;
}

KarstHole::KarstHole(const Vector3& start, const Vector3& end, float width, float height,
                     KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
    : KarstHole(width, height, startingShape, endingShape)
{
    this->path = BSpline({start, end});
}

KarstHole::KarstHole(BSpline fullPath, float width, float height,
                     KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
    : KarstHole(width, height, startingShape, endingShape)
{
    this->path = fullPath;
}

KarstHoleProfile KarstHole::interpolate(float t, float previousAcceptedTime, float nextAcceptedTime)
{
    /* Linear Interpolation, the simpliest */
    return this->startingProfile.interpolate(this->endingProfile, this->path, t, previousAcceptedTime, nextAcceptedTime); //KarstHoleProfile(interpolated);
}

std::vector<std::vector<Vector3>> KarstHole::computeClosingMesh(std::vector<Vector3>& vertices)
{
    // Compute the closing faces using the "ear clipping" method
    Vector3 middle(0, 0, 0);
    Vector3 ray; //Vector3(-1.1, -1.2, points[i].z);
    Vector3 kindOfTangent;
    std::vector<int> remaining_nodes(vertices.size()); // Starts with {0, 1, 2, 3, ...N}
    std::vector<std::vector<Vector3>> triangles;
    for (size_t i = 0; i < remaining_nodes.size(); i++) {
        middle += vertices[i];
        remaining_nodes[i] = i;
    }
    middle /= (float)number_of_points;
    for (int i = 0; i < number_of_points; i++) {
        if ((vertices[i] - middle).norm2() > ray.norm2()) {
            ray = (vertices[i] - middle) * 2.f;
            int second_indice = 1;
            while(kindOfTangent.norm2() == 0 || kindOfTangent.dot(ray) == 0) {
                kindOfTangent = (vertices[(i + second_indice)] - middle);
                second_indice++;
                if ((i + second_indice) % vertices.size() == 0) {
//                            std::cout << "Fuuuuuu..." << std::endl;
                    break;
                }
            }
        }
    }
    ray += kindOfTangent;
    ray += middle;

    int max_tries = vertices.size() * vertices.size();
    while(remaining_nodes.size() > 3) {
        max_tries --;
        if (max_tries < 0) {
            std::cout << "Breaking : too much iterations" << std::endl;
            break;
        }
        int previous = remaining_nodes[0],
                current = remaining_nodes[1],
                next = remaining_nodes[2];
        Vector3 midpoint = (vertices[previous] + vertices[next])/2.f;
        // Check if midpoint is in the (reduced) polygon
        int number_of_intersections = 0;
        bool valid = true;
        for (size_t j = 0; j < vertices.size(); j++) {
            // Count the number of intersection from the midpoint to somewhere outside
            if (Collision::intersectionBetweenTwoSegments(ray, midpoint, vertices[j], vertices[(j + 1) % vertices.size()]).isValid()) {
                number_of_intersections++;
            }
        }
        for (size_t j = 0; j < remaining_nodes.size() && valid; j++) {
            if (previous != remaining_nodes[j] && previous != remaining_nodes[(j + 1) % remaining_nodes.size()] &&
                    next != remaining_nodes[j] && next != remaining_nodes[(j + 1) % remaining_nodes.size()]) {
                // Also, check if the "previous-next" line intersects any other edge (except at the exact position of points)
                Vector3 inter = Collision::intersectionBetweenTwoSegments(vertices[previous], vertices[next], vertices[remaining_nodes[j]], vertices[remaining_nodes[(j + 1) % remaining_nodes.size()]]);
                if (inter.isValid()) {
                    if (inter != vertices[previous] && inter != vertices[next]) {
                        valid = false;

                    }
                }
            }
        }
        // If there is an odd number of itersection, point is inside the shape,
        // we can create a triangle and remove the current node
        if (valid && number_of_intersections % 2 == 1) {
            triangles.push_back({vertices[current], vertices[previous], vertices[next]});
            remaining_nodes.erase(remaining_nodes.begin() + 1);
        } else {
            // Otherwise, rotate the list just to change the first 3 values
            std::rotate(remaining_nodes.begin(), remaining_nodes.begin() + 1, remaining_nodes.end());
        }
    }
    // Put the last 3 in a triangle
    if (remaining_nodes.size() == 3) {
        triangles.push_back({vertices[remaining_nodes[1]], vertices[remaining_nodes[0]], vertices[remaining_nodes[2]]});    
    }
    return triangles;
}

std::vector<std::vector<Vector3> > KarstHole::generateMesh()
{
    this->vertexGroups.clear();
    this->vertexCylinders.clear();
    this->cylinders.clear();

    Vector3 minVec = Vector3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vector3 maxVec = Vector3(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());

    std::vector<float> v_times;
    float intervals = 1.f / (float)(int(this->path.points.size()) - 1);
    for(size_t i = 0; i < this->path.points.size() - 1; i++) {
        v_times.push_back((i +  0) * intervals);
        v_times.push_back((i + .1) * intervals);
        v_times.push_back((i + .9) * intervals);
    }
    v_times.push_back(1.f);

    number_of_intermediates = v_times.size(); //2 * int(this->path.points.size()-1) + 1; // std::min(5, int(this->path.points.size())*2);
    float dt = 1.f / (float)(number_of_intermediates - 1);
    std::vector<float> validTimesInPath;
    for (int i = 0; i < number_of_intermediates; i++) {
        float t = v_times[i]; //i * dt;
        std::vector<Vector3> intermediateShape = this->interpolate(t).vertices.points;
        for (const Vector3& pos : intermediateShape) {
            minVec.x = std::min(minVec.x, pos.x);
            minVec.y = std::min(minVec.y, pos.y);
            minVec.z = std::min(minVec.z, pos.z);
            maxVec.x = std::max(maxVec.x, pos.x);
            maxVec.y = std::max(maxVec.y, pos.y);
            maxVec.z = std::max(maxVec.z, pos.z);
        }
        // Find non-vertical times in our path
        if (!this->path.getDirection(t).isAlmostVertical()) {
            validTimesInPath.push_back(t);
        }
    }

    std::vector<std::vector<Vector3>> triangles;
    std::vector<std::vector<Vector3>> allIntermediateVertices(number_of_intermediates + 2);
    std::vector<KarstHoleProfile> intermediateProfiles;
    for (int iT = 0; iT < number_of_intermediates + 2; iT++) {
        float t = v_times[std::min(std::max(iT, 0), int(v_times.size()-1))]; //(iT-1.f) * dt;
        float previousValidTime = -1.f, nextValidTime = -1.f;
        for (const float& time : validTimesInPath) {
            if (time < t)
                previousValidTime = time;
            if (time >= t) {
                nextValidTime = time;
                break; // no need to go further
            }
        }
        KarstHoleProfile currentProfile;
        if (this->path.getDirection(t).isAlmostVertical() && iT != number_of_intermediates + 1) {
            KarstHoleProfile prev(KarstHolePredefinedShapes::TUBE);
            prev.setSize(this->width, this->height);
            prev.translate(path.getPoint(t));
//            std::cout << "Verticality at t = " << t << " as the direction is " << this->path.getDirection(t) << std::endl;
            if (iT > 0) {
//                std::cout << previousValidTime << " < " << t << " < " << nextValidTime << " -> " << interpolation::linear(t, previousValidTime, nextValidTime) << std::endl;
                currentProfile = prev; //.rotateIndicesUntilBestFitWith(intermediateProfiles.back(), this->number_of_points);
            }

        } else {
            if (iT == 0) {
                currentProfile = this->interpolate(0.0).translate(path.getDirection(0.f) * -this->width).scale(.5f);
                cylinders.push_back(std::make_tuple(path.getPoint(0.f) + path.getDirection(0.f) * -this->width, path.getPoint(0.f)));
            }
            else if(iT == number_of_intermediates + 1) {
//                std::cout << "Direction " << path.getDirection(1.f) * this->size << std::endl;
                currentProfile = this->interpolate(1.0).translate(path.getDirection(1.f) * this->width).scale(.5f);
                cylinders.push_back(std::make_tuple(path.getPoint(1.f), path.getPoint(1.f) + path.getDirection(1.f) * this->width));
            }
            else {
                currentProfile = this->interpolate(t, previousValidTime, nextValidTime);
                if (iT < number_of_intermediates) // Avoid the last cylinder
                    cylinders.push_back(std::make_tuple(path.getPoint(t), path.getPoint(t + dt)));
            }
        }
        currentProfile.vertices.closed = false;
        intermediateProfiles.push_back(currentProfile);
        allIntermediateVertices[iT] = currentProfile.vertices.getPath(number_of_points);
    }

    for (size_t i = 0; i < allIntermediateVertices.size() - 1; i++) {
        std::vector<Vector3> startingShape = intermediateProfiles[i].vertices.getPath(number_of_points);
        std::vector<Vector3> endingShape = intermediateProfiles[i + 1]/*.rotateIndicesUntilBestFitWith(intermediateProfiles[i], number_of_points)*/.vertices.getPath(number_of_points);
        // Compute the triangles that makes the actual tunnel
        for (int j = 0; j < number_of_points; j++) {
            int j_plus_1 = (j + 1) % number_of_points;
            triangles.push_back({
                                    (startingShape[j_plus_1]),
                                    (startingShape[j]),
                                    (endingShape[j])
                                });
            triangles.push_back({
                                    (endingShape[j]),
                                    (endingShape[j_plus_1]),
                                    (startingShape[j_plus_1])
                                });
        }

//        if (i == 0) {
            std::vector<std::vector<Vector3>> firstClosingShape = this->computeClosingMesh(allIntermediateVertices[i]); // For the first tunnel, compute the front shape
            triangles.insert(triangles.end(), firstClosingShape.begin(), firstClosingShape.end()); // Compute the back closing shape and add the triangles.
//        } else if (i == allIntermediateVertices.size() - 2) {
            std::vector<std::vector<Vector3>> closingShape = this->computeClosingMesh(allIntermediateVertices[i + 1]);
            triangles.insert(triangles.end(), closingShape.begin(), closingShape.end()); // Compute the back closing shape and add the triangles.
//        }

    }
    return triangles;
}

std::tuple<GridF, Vector3> KarstHole::generateMask(std::vector<std::vector<Vector3>>  precomputedTriangles)
{
    std::vector<std::vector<Vector3>> triangles = precomputedTriangles;
    if (triangles.empty())
        triangles = this->generateMesh();

    Vector3 minVec = Vector3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vector3 maxVec = Vector3(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
    for (const auto& triangle : triangles) {
        for (const Vector3& pos : triangle) {
            minVec.x = std::min(minVec.x, pos.x);
            minVec.y = std::min(minVec.y, pos.y);
            minVec.z = std::min(minVec.z, pos.z);
            maxVec.x = std::max(maxVec.x, pos.x);
            maxVec.y = std::max(maxVec.y, pos.y);
            maxVec.z = std::max(maxVec.z, pos.z);
        }
    }
    for (auto& triangle : triangles) {
        for (Vector3& pos : triangle) {
            pos -= minVec;
        }
    }
    /*
    GridF mask((maxVec - minVec) + Vector3(1, 1, 1), 0);

    for (auto& triangle : triangles) {
        Vector3 minTriangle = mask.getDimensions(), maxTriangle = Vector3(0, 0, 0);
        for (Vector3& pos : triangle) {
            minTriangle.x = std::min(minTriangle.x, pos.x);
            maxTriangle.x = std::max(maxTriangle.x, pos.x);
            minTriangle.y = std::min(minTriangle.y, pos.y);
            maxTriangle.y = std::max(maxTriangle.y, pos.y);
            minTriangle.z = std::min(minTriangle.z, pos.z);
            maxTriangle.z = std::max(maxTriangle.z, pos.z);
        }
        for (int x = minTriangle.x; x < maxTriangle.x; x++) {
            for (int y = minTriangle.y; y < maxTriangle.y; y++) {
                for (int z = minTriangle.z; z < maxTriangle.z; z++) {
                    if (Collision::intersectionTriangleAABBox(triangle[0], triangle[1], triangle[2], Vector3(x, y, z), Vector3(x+1, y+1, z+1))) {
                        mask.at(x, y, z) = 1.f;
                    }
                }
            }
        }
    }

    for (int x = 0; x < mask.sizeX; x++) {
        for (int y = 0; y < mask.sizeY; y++) {
            bool currentlyInside = false;
            bool lastVoxelWasSurface = false;
            for (int z = 0; z < mask.sizeZ; z++) {
                if (mask.at(x, y, z) > 0) {
                    if (!lastVoxelWasSurface) {
                        currentlyInside = !currentlyInside;
                    } else {
                        lastVoxelWasSurface = true;
                    }
                } else if (currentlyInside) {
                    mask.at(x, y, z) = 1.f;
                    lastVoxelWasSurface = false;
                }
            }
        }
    }

    Vector3 anchor = this->path.points[0] - minVec;
    return std::make_tuple(mask, anchor);
    */

    std::vector<Vector3> cylindersStart;
    std::vector<Vector3> cylindersEnd;
    for (auto& cylinder : this->cylinders) {
        cylindersStart.push_back(std::get<0>(cylinder) - minVec);
        cylindersEnd.push_back(std::get<1>(cylinder) - minVec);
    }
    GridF mask((maxVec - minVec));
#pragma omp parallel for collapse(3)
    for (int x = 0; x < mask.sizeX; x++) {
        for (int y = 0; y < mask.sizeY; y++) {
            for (int z = 0; z < mask.sizeZ; z++) {
                Vector3 point(x, y, z);
                bool allCollisionsValidated = false;
                int numberOfCollisions = 0;
                int currentTry = 0;
                Vector3 ray;
                while (!allCollisionsValidated) {
                    currentTry++;
                    if (currentTry > 100) {
                        numberOfCollisions = 1; // Set it inside
                        std::cout << "Way too many iterations..." << std::endl;
                        break;
                    }
                    allCollisionsValidated = true;
                    numberOfCollisions = 0;
                    ray = Vector3(-200, y, z); // + Vector3::random() * 180.f; // (Vector3::random() * 2.f * (maxVec - minVec).norm()).translate((minVec - maxVec)/ 2.f);
                    int i = 0;
                    int lastTrianglesCylinder = -1;
                    bool ignoreThisCylinder = false;
                    for (const std::vector<Vector3>& triangle : triangles) {
                        // Estimate the cylinder associated to this particular triangle.
                        int triangle_group = i / (2 * (number_of_points - 2) + 2 * number_of_points);
                        i++;
                        if (triangle_group != lastTrianglesCylinder) {
                            lastTrianglesCylinder = triangle_group;
                            if (numberOfCollisions % 2 == 1) {
                                break;
                            }
//                            ignoreThisCylinder = false;
//                            float distToCylinder = shortestDistanceBetweenSegments(point, ray, cylindersStart[triangle_group], cylindersEnd[triangle_group]);
//                            if (distToCylinder > std::max(this->width, this->height)) {
//                                ignoreThisCylinder = true;
//                            }
                        }
                        // Ignore this calculation if it's not needed.
//                        if (ignoreThisCylinder) {
//                            continue;
//                        }

                        bool collision_result = Collision::segmentToTriangleCollision(point, ray, triangle[0], triangle[1], triangle[2]).isValid();
                        if (collision_result) {
                            numberOfCollisions ++;
                        }
//                        else if (collision_result == -1) {
//                            allCollisionsValidated = false;
//                            break;
//                        }
                    }
                }
                if (numberOfCollisions % 2 == 1) {
                    mask.at(x, y, z) = 1.f;
                }
            }
        }
    }
//    mask = mask.toDistanceMap();
//    float maxDistance = mask.max();//this->size/2.f;
//    for (float& m : mask) {
//        if (m > maxDistance) m = maxDistance;
//    }
////    std::cout << "Min distances : " << mask.min() << " -- max : " << mask.max() << " -- clamped to " << maxDistance << "\n";
////    std::cout << mask.displayValues() << std::endl;
//    for (float& m : mask) {
////        if (m > 0) m = 1.f;
//        m = interpolation::linear(m, 0.f, maxDistance);
////        m = interpolation::quadratic(interpolation::linear(m, 0.f, 5.f)); //(sigmoid(m) - s_0) / (s_1 - s_0);
//    }
    Vector3 anchor = this->path.points[0] - minVec;
    return std::make_tuple(mask, anchor);
}
