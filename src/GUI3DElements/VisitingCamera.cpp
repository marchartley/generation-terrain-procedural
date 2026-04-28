#include "VisitingCamera.h"

VisitingCamera::VisitingCamera() : qglviewer::Camera()
{
    this->setFlySpeed(0);
}

void VisitingCamera::moveForward(float distance)
{
    this->setFlySpeed(0);
    std::vector<CatmullRomSpline> closestPathCandidates;
    std::vector<int> closestPointIndexCandidates;
    CatmullRomSpline closestPath;
    int closestPointIndex = -1;
    float bestScore = -1;
    Vector3 currentPos = position();
    Vector3 currentDir = viewDirection();
    for (const CatmullRomSpline& path : this->paths)
    {
        auto pathPoints = path.getPath();
        for (size_t i = 0; i < pathPoints.size(); i++)
        {
            Vector3 pointPos = pathPoints[i];
            Vector3 toPoint = pointPos - currentPos;
            if (toPoint.norm2() < 2.f) continue;

            float score = toPoint.normalized().dot(currentDir) / toPoint.norm2();
            if (score > bestScore) {
                bestScore = score;
                closestPath = path;
                closestPointIndex = i;
            }
        }
    }
    if (closestPointIndex > -1) {
        auto closestPathPoints = closestPath.getPath();
        Vector3 desiredDirection = (closestPathPoints[closestPointIndex] - currentPos).normalize();
        this->setPosition(currentPos + desiredDirection * distance);
        this->setPivotPoint(position());
//        this->lookAt(closestPath.points[closestPointIndex]);
    }
}

void VisitingCamera::moveBackward(float distance)
{
    this->setFlySpeed(0);
    CatmullRomSpline closestPath;
    float closestTime;
    float bestDist = std::numeric_limits<float>::max();
    Vector3 currentPos = position();
    Vector3 currentDir = viewDirection();
    for (CatmullRomSpline& path : this->paths)
    {
        Vector3 closestPoint = path.estimateClosestPos(currentPos);
        if ((currentPos - closestPoint).norm2() < bestDist) {
            bestDist = (currentPos - closestPoint).norm2();
            closestPath = path;
        }
    }
    if (this->paths.size() > 0) {
        closestTime = closestPath.estimateClosestTime(currentPos);
        Vector3 desiredDirection = closestPath.getDirection(closestTime) * -currentDir.dot(closestPath.getDirection(closestTime));
        this->setPosition(currentPos + desiredDirection * distance);
        this->setPivotPoint(position());
//        this->setViewDirection(desiredDirection);
    }
}
