#include "KarstHole.h"

#include "Utils/Utils.h"

KarstHole::KarstHole(float size)
{
    this->startingProfile = KarstHoleProfile(KarstHolePredefinedShapes::SOLUBLE_BED).setSize(size, size);
    this->endingProfile = KarstHoleProfile(KarstHolePredefinedShapes::KEYHOLE).setSize(size, size);
    this->path = BSpline();
    this->size = size;
}

KarstHole::KarstHole(Vector3 start, Vector3 end, float size) : KarstHole(size)
{
    this->path = BSpline({start, end});
}

KarstHole::KarstHole(BSpline fullPath, float size) : KarstHole(size)
{
    this->path = fullPath;
//    std::cout << "Start shape : \n";
//    for (const auto& pos : this->startingProfile.vertices.getPath(0.1f))
//        std::cout << "(" << pos.x << ", " << pos.y << ", " << pos.z << ") ";
//    std::cout << "\nEnding shape : \n";
//    for (const auto& pos : this->endingProfile.vertices.getPath(0.1f))
//        std::cout << "(" << pos.x << ", " << pos.y << ", " << pos.z << ") ";
}

KarstHoleProfile KarstHole::interpolate(float t)
{
    /* Linear Interpolation, the simpliest */
    return this->startingProfile.interpolate(this->endingProfile, this->path, t); //KarstHoleProfile(interpolated);
}
std::pair<KarstHoleProfile, std::vector<std::vector<Vector3>>> KarstHole::interpolateAndGetMesh(float t)
{
    /* Linear Interpolation, the simpliest */
    return this->startingProfile.interpolateAndGetMesh(this->endingProfile, this->path, t); //KarstHoleProfile(interpolated);
}

std::vector<std::vector<Vector3> > KarstHole::generateMesh()
{
    Vector3 minVec = Vector3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vector3 maxVec = Vector3(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());

    int number_of_points = 5;
    int number_of_intermediates = int(this->path.points.size())*2; // std::min(5, int(this->path.points.size())*2);
    float dt = 1.f / (float)(number_of_intermediates - 1);
    for (int i = 0; i < number_of_intermediates; i++) {
        float t = i * dt;
        std::vector<Vector3> intermediateShape = this->interpolate(t).vertices.points;
        for (const Vector3& pos : intermediateShape) {
            minVec.x = std::min(minVec.x, pos.x);
            minVec.y = std::min(minVec.y, pos.y);
            minVec.z = std::min(minVec.z, pos.z);
            maxVec.x = std::max(maxVec.x, pos.x);
            maxVec.y = std::max(maxVec.y, pos.y);
            maxVec.z = std::max(maxVec.z, pos.z);
        }
    }

    std::vector<std::vector<Vector3>> triangles;
    std::vector<std::vector<Vector3>> allIntermediateVertices(number_of_intermediates + 2);
    for (int iT = 0; iT < number_of_intermediates + 2; iT++) {
        float t = (iT-1.f) * dt;
        std::vector<Vector3> points;
        if (iT == 0)
            points = this->interpolate(0.0).translate(path.getDirection(0.f) * -this->size).scale(.5f).vertices.getPath(1.f / (float)(number_of_points - 1));
        else if(iT == number_of_intermediates + 1)
            points = this->interpolate(1.0).translate(path.getDirection(1.f) * this->size).scale(.5f).vertices.getPath(1.f / (float)(number_of_points - 1));
        else
            points = this->interpolate(t).vertices.getPath(1.f / (float)(number_of_points - 1));

        allIntermediateVertices[iT] = points;
        // Add the closing faces only on the first and last of the path using the "ear clipping" method
        if (iT == 0 || iT == number_of_intermediates + 1) {
            Vector3 middle(0, 0, 0);
            Vector3 ray; //Vector3(-1.1, -1.2, points[i].z);
            Vector3 kindOfTangent;
            std::vector<int> remaining_nodes(number_of_points); // Starts with {0, 1, 2, 3, ...N}
            for (int i = 0; i < number_of_points; i++) {
                middle += points[i];
                remaining_nodes[i] = i;
            }
            middle /= (float)number_of_points;
            for (int i = 0; i < number_of_points; i++) {
                if ((points[i] - middle).norm2() > ray.norm2()) {
                    ray = (points[i] - middle) * 2.f;
                    int second_indice = 1;
                    while(kindOfTangent.norm2() == 0 || kindOfTangent.dot(ray) == 0) {
                        kindOfTangent = (points[(i + second_indice)] - middle);
                        second_indice++;
                        if ((i + second_indice) % points.size() == 0) {
//                            std::cout << "Fuuuuuu..." << std::endl;
                            break;
                        }
                    }
                }
            }
            ray += kindOfTangent;
            ray += middle;

            int max_tries = number_of_points * number_of_points;
            while(remaining_nodes.size() > 3) {
                max_tries --;
                if (max_tries < 0) {
                    std::cout << "Breaking : too much iterations" << std::endl;
                    for (int i = 0; i < number_of_points; i++) {
                        std::cout << points[i] << "\n";
                    }
                    std::cout << std::endl;
                    break;
                }
                int previous = remaining_nodes[0],
                        current = remaining_nodes[1],
                        next = remaining_nodes[2];
                Vector3 midpoint = (points[previous] + points[next])/2.f;
                // Check if midpoint is in the (reduced) polygon
                int number_of_intersections = 0;
                bool valid = true;
                int intersectionId1, intersectionId2;
                for (size_t j = 0; j < points.size(); j++) {
                    // Count the number of intersection from the midpoint to somewhere outside
                    if (intersection(ray, midpoint, points[j], points[(j + 1) % points.size()])) {
                        number_of_intersections++;
                    }
                }
                for (size_t j = 0; j < remaining_nodes.size() && valid; j++) {
                    if (previous != remaining_nodes[j] && previous != remaining_nodes[(j + 1) % remaining_nodes.size()] &&
                            next != remaining_nodes[j] && next != remaining_nodes[(j + 1) % remaining_nodes.size()]) {
                        // Also, check if the "previous-next" line intersects any other edge (except at the exact position of points)
                        if (intersection(points[previous], points[next], points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]])) {
                            Vector3 inter = intersectionPoint(points[previous], points[next], points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]]);
                            if (inter != points[previous] && inter != points[next]) {
                                valid = false;
                                intersectionId1 = remaining_nodes[j];
                                intersectionId2 = remaining_nodes[(j + 1) % remaining_nodes.size()];

                            }
                        }
                    }
                }
                // If there is an odd number of itersection, point is inside the shape,
                // we can create a triangle and remove the current node
                if (valid && number_of_intersections % 2 == 1) {
                    triangles.push_back({points[previous], points[current], points[next]});
                    remaining_nodes.erase(remaining_nodes.begin() + 1);
                } else {
                    // Otherwise, rotate the list just to change the first 3 values
                    std::rotate(remaining_nodes.begin(), remaining_nodes.begin() + 1, remaining_nodes.end());
                    if (max_tries < 10) {
                        std::cout << "For polygon " << iT << " and triangle " << previous << "-" << current << "-" << next << ": ";
                        if (!valid)
                            std::cout << "\n- The line " << previous << "-" << next << " has an intersection with line " << intersectionId1 << "-" << intersectionId2 << ".";
                        if (number_of_intersections % 2 == 0)
                            std::cout << "\n- The midpoint (" << midpoint << ") is outside the shape (" << number_of_intersections << " intersections)";
                        std::cout << "\n--> middle is at (" << middle << ") and ray is " << ray << std::endl;
//                        std::cout << std::endl;
                    }
                }
            }
            // Put the last 3 in a triangle
            if (remaining_nodes.size() == 3)
                triangles.push_back({points[remaining_nodes[0]], points[remaining_nodes[1]], points[remaining_nodes[2]]});
        }
    }

    for (int i = 0; i < number_of_intermediates + 1; i++) {
        for (int j = 0; j < number_of_points; j++) {
            int j_plus_1 = (j + 1) % number_of_points;
            triangles.push_back({
                                    (allIntermediateVertices[i][j]),
                                    (allIntermediateVertices[i][j_plus_1]),
                                    (allIntermediateVertices[i+1][j])
                                });
            triangles.push_back({
                                    (allIntermediateVertices[i+1][j]),
                                    (allIntermediateVertices[i+1][j_plus_1]),
                                    (allIntermediateVertices[i][j_plus_1])
                                });
        }
    }
    return triangles;
}

std::tuple<Matrix3<float>, Vector3> KarstHole::generateMask(std::vector<std::vector<Vector3>>  precomputedTriangles)
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
    Matrix3<float> mask((maxVec - minVec)); //Vector3((maxVec - minVec).x, (maxVec - minVec).y, 1));
    for (int x = 0; x < mask.sizeX; x++) {
        for (int y = 0; y < mask.sizeY; y++) {
            for (int z = 0; z < mask.sizeZ; z++) {
                Vector3 point(x, y, z);
                bool allCollisionsValidated = false;
                int numberOfCollisions = 0;
                int currentTry = 0;
                while (!allCollisionsValidated) {
                    currentTry++;
                    if (currentTry > 100) {
                        numberOfCollisions = 1; // Set it inside
                        std::cout << "Way too much iterations..." << std::endl;
                        break;
                    }
                    allCollisionsValidated = true;
                    numberOfCollisions = 0;
                    Vector3 ray = (Vector3::random().normalize() * 200.f * maxVec).translate((maxVec + minVec)/ 2.f);
                    for (const std::vector<Vector3>& triangle : triangles) {
                        int collision_result = this->segmentToTriangleCollision(point, ray, triangle[0], triangle[1], triangle[2]);
                        if (collision_result == 1)
                            numberOfCollisions ++;
                        else if (collision_result == -1) {
                            allCollisionsValidated = false;
                            break;
                        }
                    }
                }
                if (numberOfCollisions % 2 == 1) {
                    mask.at(x, y, z) = 1.f;
                }
            }
        }
    }
    mask = mask.toDistanceMap();
    float maxDistance = mask.max();//this->size/2.f;
    /*for (float& m : mask) {
        if (m > maxDistance) m = maxDistance;
    }*/
//    std::cout << "Min distances : " << mask.min() << " -- max : " << mask.max() << " -- clamped to " << maxDistance << "\n";
//    std::cout << mask.displayValues() << std::endl;
    for (float& m : mask) {
        m = interpolation::linear(m, 0.f, maxDistance);
//        m = interpolation::quadratic(interpolation::linear(m, 0.f, 5.f)); //(sigmoid(m) - s_0) / (s_1 - s_0);
    }
    Vector3 anchor = this->path.points[0] - minVec;
    return std::make_tuple(mask, anchor);
}

float tetrahedronSignedVolume(Vector3 a, Vector3 b, Vector3 c, Vector3 d) {
    return (1.f/6.f) * (c-a).cross((b-a)).dot((d-a));
}
int sign(float a){return (a < 0 ? -1 : (a > 0 ? 1 : 0)); }
int KarstHole::segmentToTriangleCollision(Vector3 s1, Vector3 s2, Vector3 t1, Vector3 t2, Vector3 t3)
{
    Vector3 ray = (s1 - s2).normalize();
    Vector3 segment1 = (s1 - t1).normalize();
    Vector3 segment2 = (s1 - t2).normalize();
    Vector3 segment3 = (s1 - t3).normalize();
    if (ray.dot(segment1) > 0.99 || ray.dot(segment2) > 0.99 || ray.dot(segment3) > 0.99) {
//        std::cout << ray << "." << segment1 << " = " << ray.dot(segment1) << " " << ray << "." << segment2 << " = " << ray.dot(segment2) << " " << ray << "." << segment3 << " = " << ray.dot(segment3) << std::endl;
        return -1;
    }
    // Compute tetraheadron "signed" volume from one end of the segment and the triangle
    float product_of_volumes = tetrahedronSignedVolume(s1, t1, t2, t3) * tetrahedronSignedVolume(s2, t1, t2, t3);
    if (std::abs(product_of_volumes) < 0.000001 || sign(product_of_volumes) == 0)
        return -1;
    if (sign(product_of_volumes) > 0) // Same sign for the two volums computed
        return 0;
    // Compute the tetrahedron built with the segment
    // and 2 points of the triangle, all configs should have
    // the same sign.
    int volume1 = tetrahedronSignedVolume(s1, s2, t1, t2);
    int volume2 = tetrahedronSignedVolume(s1, s2, t3, t1);
    int volume3 = tetrahedronSignedVolume(s1, s2, t2, t3);
    if (std::abs(volume1) < 0.000001 || std::abs(volume2) < 0.000001 || std::abs(volume3) < 0.000001)
        return -1;
    if (sign(volume1) == sign(volume2) && sign(volume1) == sign(volume3))
        return 1;
    return 0;
}
