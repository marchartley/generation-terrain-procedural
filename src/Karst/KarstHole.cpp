#include "KarstHole.h"

#include "Utils/Utils.h"

KarstHole::KarstHole()
{
    this->startingProfile = KarstHoleProfile(KarstHolePredefinedShapes::SOLUBLE_BED).setSize(15, 15);
    this->endingProfile = KarstHoleProfile(KarstHolePredefinedShapes::KEYHOLE).setSize(15, 15);
    this->path = BSpline();
}

KarstHole::KarstHole(Vector3 start, Vector3 end) : KarstHole()
{
    this->path = BSpline({start, end});
}

KarstHole::KarstHole(BSpline fullPath) : KarstHole()
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
    std::vector<Vector3> startingPoints = this->startingProfile.vertices.getPath(.1f);
    std::vector<Vector3> endingPoints = this->endingProfile.vertices.getPath(.1f);
    std::vector<Vector3> interpolated(startingPoints.size());

    for (size_t i = 0; i < startingPoints.size(); i++) {
        interpolated[i] = (startingPoints[i] * (1 - t) + endingPoints[i] * t).setDirection(this->path.getDerivative(t)).translate(this->path.getPoint(t));
    }
    return KarstHoleProfile(interpolated);
}

std::tuple<Matrix3<float>, Vector3> KarstHole::generateMask()
{
    Vector3 minVec = Vector3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vector3 maxVec = Vector3(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());

    int number_of_points = 10;
    int number_of_intermediates = 5;
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

    Matrix3<float> mask((maxVec - minVec)); //Vector3((maxVec - minVec).x, (maxVec - minVec).y, 1));

    std::vector<std::vector<Vector3>> triangles;
    std::vector<std::vector<Vector3>> allIntermediateVertices(number_of_intermediates);
    for (int i = 0; i < number_of_intermediates; i++) {
        float t = i * dt;
        std::vector<Vector3> points = this->interpolate(t).vertices.getPath(1.f / (float)(number_of_points - 1));
        for (Vector3& p : points) {
            p -= minVec;
            std::cout << p << " ";
        }
        std::cout << std::endl;
        allIntermediateVertices[i] = points;
        // Add the closing faces only on the first and last of the path using the "ear clipping" method
        if (i == 0 || i == number_of_intermediates - 1) {
            Vector3 ray;
//            if (i == 0)
//            ray = Vector3(minVec.x, minVec.y, points[i].z);
            ray = Vector3(-1.1, -1.2, points[i].z);
            std::vector<int> remaining_nodes(number_of_points); // Starts with {0, 1, 2, 3, ...N}
            for (int i = 0; i < number_of_points; i++) remaining_nodes[i] = i;
            int max_tries = number_of_points * number_of_points;
            while(remaining_nodes.size() > 3) {
                max_tries --;
                if (max_tries < 0) break;
//                for (size_t j = 0; j < remaining_nodes.size(); j++)
//                    std::cout << remaining_nodes[j] << " ";
                int previous = remaining_nodes[0],
                        current = remaining_nodes[1],
                        next = remaining_nodes[2];
                Vector3 midpoint = (points[previous] + points[next])/2.f;
//                std::cout << midpoint << std::endl;
                // Check if midpoint is in the (reduced) polygon
                int number_of_intersections = 0;
                bool valid = true;
                for (size_t j = 0; j < remaining_nodes.size() && valid; j++) {
                    // Count the number of intersection from the midpoint to somewhere outside
                    if (intersection(ray, midpoint, points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]])) {
                        number_of_intersections++;
//                        std::cout << "Intersection between nodes (" << minVec << ", " << midpoint << ") and (" << remaining_nodes[j] << ", " << remaining_nodes[(j + 1) % remaining_nodes.size()] << ")\n";
                    } else {
//                        std::cout << "No intersection between nodes (" << minVec << ", " << midpoint << ") and (" << remaining_nodes[j] << ", " << remaining_nodes[(j + 1) % remaining_nodes.size()] << ")\n";
                    }
                    // Also, check if the "previous-next" line intersects any other edge (except at the exact position of points)
                    if (intersection(points[previous], points[next], points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]])) {
                        Vector3 inter = intersectionPoint(points[previous], points[next], points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]]);
                        if (inter != points[previous] && inter != points[next]) {
                            valid = false;
//                            std::cout << "Intersection between nodes (" << previous << ", " << next << ") and (" << remaining_nodes[j] << ", " << remaining_nodes[(j + 1) % remaining_nodes.size()] << ")\n";
                        } else {
//                            std::cout << "Exactly on the nodes (" << previous << ", " << next << ") and (" << remaining_nodes[j] << ", " << remaining_nodes[(j + 1) % remaining_nodes.size()] << ")\n";
                        }
                    } else {
//                        std::cout << "Didn't even touch nodes (" << previous << ", " << next << ") and (" << remaining_nodes[j] << ", " << remaining_nodes[(j + 1) % remaining_nodes.size()] << ")\n";
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
                }
            }
            // Put the last 3 in a triangle
            if (remaining_nodes.size() == 3)
                triangles.push_back({points[remaining_nodes[0]], points[remaining_nodes[1]], points[remaining_nodes[2]]});
        }
    }
//    std::cout << "We have " << triangles.size() << " triangles and we should have " << 2*(number_of_points-2) << std::endl;
    for (int i = 0; i < number_of_intermediates - 1; i++) {
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
    std::cout << "Triangles : \n";
    for (const auto& t : triangles) {
        std::cout << t[0] << " " << t[1] << " " << t[2] << "\n";
    }

    for (int x = 0; x < mask.sizeX; x++) {
        for (int y = 0; y < mask.sizeY; y++) {
            for (int z = 0; z < mask.sizeZ; z++) {
                Vector3 point(x, y, z);
                Vector3 ray(-110, -120, -130);
                int numberOfCollisions = 0;
                for (const std::vector<Vector3>& triangle : triangles) {
                    if (this->segmentToTriangleCollision(point, ray, triangle[0], triangle[1], triangle[2]))
                        numberOfCollisions ++;
//                    std::cout << "(" << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << ") ";
                }
                if (numberOfCollisions % 2 == 1) {
                    mask.at(x, y, z) = 1.f;
                }
//                std::cout << numberOfCollisions << " ";
            }
        }
    }
//    std::cout << "\n" << mask;
//    std::cout << std::endl;
    Vector3 anchor = this->path.points[0] - minVec;
    return std::make_tuple(mask, anchor);
}

float tetrahedronSignedVolume(Vector3 a, Vector3 b, Vector3 c, Vector3 d) {
    return (1.f/6.f) * (c-a).cross((b-a)).dot((d-a));
}
int sign(float a){return (a < 0 ? -1 : (a > 0 ? 1 : 0)); }
bool KarstHole::segmentToTriangleCollision(Vector3 s1, Vector3 s2, Vector3 t1, Vector3 t2, Vector3 t3)
{
    // Compute tetraheadron "signed" volume from one end of the segment and the triangle
    bool segmentsEndOnDifferentSides = sign(tetrahedronSignedVolume(s1, t1, t2, t3)) != sign(tetrahedronSignedVolume(s2, t1, t2, t3));
    if (!segmentsEndOnDifferentSides)
        return false;
    // Compute the tetrahedron built with the segment
    // and 2 points of the triangle, all configs should have
    // the same sign.
    int volume1 = sign(tetrahedronSignedVolume(s1, s2, t1, t2));
    int volume2 = sign(tetrahedronSignedVolume(s1, s2, t3, t1));
    int volume3 = sign(tetrahedronSignedVolume(s1, s2, t2, t3));
    if (volume1 == volume2 && volume1 == volume3)
        return true;
    return false;
}
