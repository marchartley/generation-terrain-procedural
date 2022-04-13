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
    case CRACK:
        this->vertices = KarstHoleProfile::createCrackProfile();
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
/*
    Vector3 new_dir = path.getDirection(t);
    Vector3 forward(0, 1, 0);
    Vector3 up(0, 0, 1);
    Vector3 right(1, 0, 0);
//    right = path.getBinormal(t);
//    up = path.getNormal(t);

//    float angle = std::acos(forward.dot(new_dir));
//    Vector3 rotation = new_dir.cross(forward).normalize();
//    right.rotate(angle, rotation);
    if (std::abs(new_dir.dot(Vector3(0, 0, 1))) < 0.999)
        right = Vector3(0, 0, 1).cross(new_dir).normalize();
    else
        right = Vector3(0, 1, 1).cross(new_dir).normalize();
    up = new_dir.cross(right).normalize();
//    std::cout << new_dir << " " << right << " " << up << std::endl;

    std::cout << "t = " << t << "\n";
    std::cout << right << " against " << path.getFrenetBinormal(t) << "\n";
    std::cout << new_dir << " against " << path.getFrenetDirection(t) << "\n";
    std::cout << up << " against " << path.getFrenetNormal(t) << std::endl;*/
    for (Vector3& point : this->vertices.points)
        point.changeBasis(path.getFrenetBinormal(t), path.getFrenetDirection(t), path.getFrenetNormal(t));
//        point.changeBasis(right, new_dir, up);
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

KarstHoleProfile KarstHoleProfile::interpolate(KarstHoleProfile other, BSpline path, float t, float previousAcceptedTime, float nextAcceptedTime)
{
    /*std::vector<Vector3> startingPoints = this->vertices.getPath(.1f);
    std::vector<Vector3> endingPoints = other.vertices.getPath(.1f);
    std::vector<Vector3> interpolatedPoints(startingPoints.size());
    for (size_t i = 0; i < startingPoints.size(); i++) {
        interpolatedPoints[i] = (startingPoints[i] * (1 - t) + endingPoints[i] * t);
    }
    KarstHoleProfile interpolation(interpolatedPoints);
    return interpolation.rotateTowardVector(path, t).translate(path.getPoint(t));
    */
    if (!path.getDirection(t).isAlmostVertical()) {
        // Non-vertical path
        std::vector<Vector3> startingPoints = this->vertices.getPath(.1f);
        std::vector<Vector3> endingPoints = other.vertices.getPath(.1f);
        std::vector<Vector3> interpolatedPoints(startingPoints.size());
        for (size_t i = 0; i < startingPoints.size(); i++) {
            interpolatedPoints[i] = (startingPoints[i] * (1 - t) + endingPoints[i] * t);
        }
        KarstHoleProfile interpolation(interpolatedPoints);
        return interpolation.rotateTowardVector(path, t).translate(path.getPoint(t));
    } else {
        // Vertical case
        KarstHoleProfile interpolation(KarstHolePredefinedShapes::TUBE);
        interpolation.setSize(this->scaling.x, this->scaling.y);
        //return interpolation.rotateTowardVector(path, t).translate(path.getPoint(t));
        if (previousAcceptedTime > -1.f && nextAcceptedTime > -1.f) {
            float binormalRotation = 0.f;
            Vector3 prevBinormal = path.getFrenetBinormal(previousAcceptedTime);
            Vector3 nextBinormal = path.getFrenetBinormal(nextAcceptedTime);
            Vector3 prevDir = path.getFrenetDirection(previousAcceptedTime);
            Vector3 nextDir = path.getFrenetDirection(nextAcceptedTime);
            Vector3 prevNorm = path.getFrenetNormal(previousAcceptedTime);
            Vector3 nextNorm = path.getFrenetNormal(nextAcceptedTime);
            Vector3 rotator = nextDir.cross(prevDir);
            binormalRotation = std::acos(prevDir.dot(nextDir));
//            binormalRotation = std::acos(nextNorm.dot(prevBinormal));
            binormalRotation *= interpolation::linear(t, previousAcceptedTime, nextAcceptedTime);
            std::cout << "t-1 = " << previousAcceptedTime << " t+1 = " << nextAcceptedTime << " -> dot = " << prevBinormal.dot(nextBinormal) << "\n";
            std::cout << "-> normal dot : " << prevNorm.dot(nextNorm) << std::endl;
            // Get the same as the previous timestamp
//            interpolation.rotateTowardVector(path, previousAcceptedTime);
            // Increase the rotation to fit the interpolation
            for (auto& point : interpolation.vertices.points) {
                if (prevBinormal.dot(nextBinormal) < 0){
                    point.rotate(3.141592, 0, 0);
                }

                point.rotate(binormalRotation, rotator);
//                point.rotate(0, 0, binormalRotation);
            }
        }
        // The points are automatically vertically-defined
        interpolation.setSize(this->scaling.x, this->scaling.y);
        return interpolation.translate((path.getPoint(t)));
    }
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

KarstHoleProfile &KarstHoleProfile::rotateIndicesUntilBestFitWith(KarstHoleProfile &otherProfile, int numberOfPointsUsed)
{
    std::vector<Vector3> startPositions = this->vertices.getPath(1.f / (float) numberOfPointsUsed);
    std::vector<Vector3> endPositions = otherProfile.vertices.getPath(1.f / (float) numberOfPointsUsed);
    float minDist = -1;
    int bestRotateFit = -1;
    bool best_reversed = false;

    for (int rot = 0; rot < numberOfPointsUsed; rot++) {
        float currentDist = 0;
        bool need_to_be_reversed = false;
        for (int i_init = 0; i_init < numberOfPointsUsed; i_init++) {
            int i = (i_init + rot) % numberOfPointsUsed;
            float t = (i / (float)(numberOfPointsUsed - 1));
            float t_init = (i_init / (float)(numberOfPointsUsed - 1));
            Vector3 moveVec = startPositions[i] - endPositions[i_init];
            currentDist += moveVec.norm2();
        }
        if (currentDist < minDist || bestRotateFit == -1) {
            minDist = currentDist;
            bestRotateFit = rot;
            best_reversed = need_to_be_reversed;
        }
    }

    for (int rot = 0; rot < numberOfPointsUsed; rot++) {
        float currentDist = 0;
        bool need_to_be_reversed = true;
        for (int i_init = 0; i_init < numberOfPointsUsed; i_init++) {
            int i = ((numberOfPointsUsed - i_init - 1) + rot) % numberOfPointsUsed;
            float t = (i / (float)(numberOfPointsUsed - 1));
            float t_init = (i_init / (float)(numberOfPointsUsed - 1));
            Vector3 moveVec = startPositions[i] - endPositions[i_init];
            currentDist += moveVec.norm2();
        }
        if (currentDist < minDist || bestRotateFit == -1) {
            minDist = currentDist;
            bestRotateFit = rot;
            best_reversed = need_to_be_reversed;
        }
    }

    // bestRotateFit = src_shape.length + bestRotateFit * (best_reversed ? 1 : -1);
    std::vector<Vector3> newPoints;
    for (int i = bestRotateFit; i < numberOfPointsUsed + bestRotateFit; i++) {
        newPoints.insert((best_reversed ? newPoints.begin() : newPoints.end() ), startPositions[i % numberOfPointsUsed]);
    }
    this->vertices.points = newPoints;
    return *this;
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
    Vector3 previousScaling = this->scaling;
    scaling = Vector3(sizeX, sizeY, 1.f); // Dunno what to do...
    for (Vector3& point : this->vertices.points)
        point *= (scaling / previousScaling);
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
            if (intersectionBetweenTwoSegments(ray, midpoint, points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]])) {
                number_of_intersections++;
            }
            // Also, check if the "previous-next" line intersects any other edge (except at the exact position of points)
            Vector3 inter = intersectionBetweenTwoSegments(points[previous], points[next], points[remaining_nodes[j]], points[remaining_nodes[(j + 1) % remaining_nodes.size()]]);
            if (inter.isValid()) {
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
                   })/*.close()*/;
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
                   })/*.close()*/;
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
                   })/*.close()*/;
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
                   })/*.close()*/;
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
                   })/*.close()*/;
}
BSpline KarstHoleProfile::createCrackProfile()
{
    return BSpline({
                       {-.25,  0, 0},
                       {-.2,  .1, 0},
                       { .2,  .1, 0},
                       { .25,  0, 0},
                       {  0, -1., 0}
                   })/*.close()*/;
}
