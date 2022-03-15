#include "KarstHoleProfile.h"
#include "Utils/Utils.h"

KarstHoleProfile::KarstHoleProfile()
{
    this->vertices = KarstHoleProfile::createTubeProfile();
}

KarstHoleProfile::KarstHoleProfile(KarstHolePredefinedShapes shape)
{
    switch (shape) {
    case TUBE:
        this->vertices = KarstHoleProfile::createTubeProfile();
        break;
    case SOLUBLE_BED:
        this->vertices = KarstHoleProfile::createSolubleBedProfile();
        break;
    case PASSAGE:
        this->vertices = KarstHoleProfile::createPassageProfile();
        break;
    case KEYHOLE:
        this->vertices = KarstHoleProfile::createKeyholeProfile();
        break;
    case CANYON:
        this->vertices = KarstHoleProfile::createCanyonProfile();
        break;
    }
}

KarstHoleProfile::KarstHoleProfile(BSpline shape)
{
    this->vertices = shape;
}

KarstHoleProfile::KarstHoleProfile(std::vector<Vector3> shape)
{
    this->vertices = BSpline(shape);
}

KarstHoleProfile &KarstHoleProfile::rotateTowardVector(BSpline path, float t)
{
    for (Vector3& point : this->vertices.points)
        point.rotate(-3.141592/2.f, 0, 0);
    Vector3 new_dir = path.getDirection(t);
    Vector3 forward(0, 1, 0);
    Vector3 up(0, 0, 1);
    Vector3 right(1, 0, 0);
    right = path.getBinormal(t);
    up = path.getNormal(t);
    /*
    float angle = std::acos(forward.dot(new_dir));
    Vector3 rotation = new_dir.cross(forward).normalize();
    right.rotate(angle, rotation);
    if (std::abs(new_dir.dot(Vector3(0, 0, 1))) < 0.999)
        right = Vector3(0, 0, 1).cross(new_dir);
    else
        right = Vector3(0, 1, 1).cross(new_dir).normalize();
    up = new_dir.cross(right);
//    std::cout << new_dir << " " << right << " " << up << std::endl;
*/
    for (Vector3& point : this->vertices.points)
        point.changeBasis(right, new_dir, up);
    return *this;
}

KarstHoleProfile &KarstHoleProfile::translate(Vector3 translation, bool verbose)
{
    if (verbose)
        std::cout << "Points are going to move " << translation << "\n";
    for (Vector3& point : this->vertices.points) {
        if (verbose)
            std::cout << " -> " << point << " -> ";
        point.translate(translation);
        if (verbose)
            std::cout << point << "\n";
    }
    if (verbose)
        std::cout << std::endl;
    return *this;
}

KarstHoleProfile &KarstHoleProfile::scale(float scale_factor, bool verbose)
{
    Vector3 mean;
    for (const Vector3& point : this->vertices.points)
        mean += point;
    if (verbose)
        std::cout << "For scale, mean = " << mean << "/" << (float)this->vertices.points.size() << " = ";
    mean /= (float)this->vertices.points.size();
    if (verbose)
        std::cout << mean << "\nPoints are going towards mean.\n";
    for (Vector3& point : this->vertices.points) {
        if (verbose)
            std::cout << " -> point " << point << " moves toward " << (mean - point) << "*" << scale_factor;
        point += (mean - point) * scale_factor;
        if (verbose)
            std::cout << " ==> " << point << "\n";
    }
    return *this;
}

KarstHoleProfile KarstHoleProfile::interpolate(KarstHoleProfile other, BSpline path, float t)
{
    std::vector<Vector3> startingPoints = this->vertices.getPath(.1f);
    std::vector<Vector3> endingPoints = other.vertices.getPath(.1f);
    std::vector<Vector3> interpolatedPoints(startingPoints.size());
    for (size_t i = 0; i < startingPoints.size(); i++) {
        interpolatedPoints[i] = (startingPoints[i] * (1 - t) + endingPoints[i] * t);
    }
    KarstHoleProfile interpolation(interpolatedPoints);
//    std::cout << "Path length : " << path.length() << " and t = " << t << " -> derivative : " << path.getDerivative(t) << "\n";
    return interpolation.rotateTowardVector(path, t).translate(path.getPoint(t));
}
std::pair<KarstHoleProfile, std::vector<std::vector<Vector3>>> KarstHoleProfile::interpolateAndGetMesh(KarstHoleProfile other, BSpline path, float t)
{
    std::vector<Vector3> startingPoints = this->vertices.getPath(.1f);
    std::vector<Vector3> endingPoints = other.vertices.getPath(.1f);
    std::vector<Vector3> interpolatedPoints(startingPoints.size());
    for (size_t i = 0; i < startingPoints.size(); i++) {
        interpolatedPoints[i] = (startingPoints[i] * (1 - t) + endingPoints[i] * t);
    }
    KarstHoleProfile interpolation(interpolatedPoints);
    std::vector<std::vector<int>> triangleIndices = this->computeTrianglesIndices(interpolatedPoints);
    interpolation = interpolation.rotateTowardVector(path, t).translate(path.getPoint(t));
    std::vector<std::vector<Vector3>> triangles;
    for (const auto& triangle : triangleIndices) {
        triangles.push_back({interpolatedPoints[triangle[0]], interpolatedPoints[triangle[1]], interpolatedPoints[triangle[2]]});
    }
    return std::make_pair(interpolation, triangles);
}

KarstHoleProfile& KarstHoleProfile::setNumberOfVertices(int vertice_count)
{
    std::vector<Vector3> newPoints = this->vertices.getPath(1/(float)(vertice_count - 1));
    this->vertices = BSpline(newPoints);
    return *this;
}

KarstHoleProfile &KarstHoleProfile::setSize(float sizeX, float sizeY)
{
    /*
    Vector3 minVec(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 0);
    Vector3 maxVec(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), 1);
    for (const Vector3& point : this->vertices.points) {
        minVec.x = std::min(minVec.x, point.x);
        minVec.y = std::min(minVec.y, point.y);
        maxVec.x = std::max(maxVec.x, point.x);
        maxVec.y = std::max(maxVec.y, point.y);
    }
    Vector3 scaling = Vector3(sizeX, sizeY, 1) / (maxVec - minVec);*/
    float scaling = std::max(sizeX, sizeY); // Dunno what to do...
    for (Vector3& point : this->vertices.points)
        point *= scaling;
    return *this;
}

std::vector<std::vector<int> > KarstHoleProfile::computeTrianglesIndices(const std::vector<Vector3>& points)
{
    std::vector<std::vector<int> > triangles;
    /// For now, only compute them when not transformed (rotated/translated)
    int number_of_points = points.size();
    Vector3 ray = Vector3(-1.1, -1.2, 0);
    std::vector<int> remaining_nodes(number_of_points); // Starts with {0, 1, 2, 3, ...N}
    for (int i = 0; i < number_of_points; i++) remaining_nodes[i] = i;
    int max_tries = number_of_points * number_of_points;
    while(remaining_nodes.size() > 3) {
        max_tries --;
        if (max_tries < 0) {
            std::cout << "Breaking : too much iterations" << std::endl;
            break;
        }
        int previous = remaining_nodes[0],
                current = remaining_nodes[1],
                next = remaining_nodes[2];
        Vector3 midpoint = (points[previous] + points[next])/2.f;
        // Check if midpoint is in the (reduced) polygon
        int number_of_intersections = 0;
        bool valid = true;
        for (size_t j = 0; j < remaining_nodes.size() && valid; j++) {
            // Count the number of intersection from the midpoint to somewhere outside
            if (intersection(ray, midpoint, points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]])) {
                number_of_intersections++;
            }
            // Also, check if the "previous-next" line intersects any other edge (except at the exact position of points)
            if (intersection(points[previous], points[next], points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]])) {
                Vector3 inter = intersectionPoint(points[previous], points[next], points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]]);
                if (inter != points[previous] && inter != points[next]) {
                    valid = false;
                }
            }
        }
        // If there is an odd number of itersection, point is inside the shape,
        // we can create a triangle and remove the current node
        if (valid && number_of_intersections % 2 == 1) {
            triangles.push_back({previous, current, next});
            remaining_nodes.erase(remaining_nodes.begin() + 1);
        } else {
            // Otherwise, rotate the list just to change the first 3 values
            std::rotate(remaining_nodes.begin(), remaining_nodes.begin() + 1, remaining_nodes.end());
            if (max_tries < 10) {
                std::cout << "For triangle " << previous << "-" << current << "-" << next << ": ";
                if (!valid)
                    std::cout << "\n- The line " << previous << "-" << next << " has an intersection.";
                if (number_of_intersections % 2 == 0)
                    std::cout << "\n- The midpoint (" << midpoint << ") is outside the shape (" << number_of_intersections << " intersections)";
                std::cout << std::endl;
            }
        }
    }
    // Put the last 3 in a triangle
    if (remaining_nodes.size() == 3)
        triangles.push_back({remaining_nodes[0], remaining_nodes[1], remaining_nodes[2]});
    return triangles;
}









BSpline KarstHoleProfile::createTubeProfile()
{
    return BSpline({
                       {-1.,  .0, 0},
                       {-.7,  .7, 0},
                       { .0,  1., 0},
                       { .7,  .7, 0},
                       { 1.,  .0, 0},
                       { .7, -.7, 0},
                       { .0, -1., 0},
                       {-.7, -.7, 0}
                   }); //.close();
}

BSpline KarstHoleProfile::createSolubleBedProfile()
{
    return BSpline({
                       {-1.,  0., 0},
                       {-.7,  .5, 0},
                       { .7,  .5, 0},
                       { 1.,  0., 0},
                       { .7, -.5, 0},
                       {-.7, -.5, 0}
                   }); //.close();
}

BSpline KarstHoleProfile::createPassageProfile()
{
    return BSpline({
                       {-1.,  0., 0},
                       {-.5,  .5, 0},
                       { .5,  .5, 0},
                       { 1.,  0., 0},
                       { .5, -.5, 0},
                       {-.5, -.5, 0}
                   }); //.close();
}

BSpline KarstHoleProfile::createKeyholeProfile()
{
    return BSpline({
                       {-1.,  0., 0},
                       {-.5,  .5, 0},
                       { .5,  .5, 0},
                       { 1.,  0., 0},
                       { .5, -.5, 0},
                       { .2, -1., 0},
                       {-.2, -1., 0},
                       {-.5, -.5, 0}
                   }); //.close();
}

BSpline KarstHoleProfile::createCanyonProfile()
{
    return BSpline({
                       {-.5,  .5, 0},
                       {-.2,  1., 0},
                       { .2,  1., 0},
                       { .5,  .5, 0},
                       { .5, -.5, 0},
                       { .2, -1., 0},
                       {-.2, -1., 0},
                       {-.5, -.5, 0}
                   }); //.close();
}
