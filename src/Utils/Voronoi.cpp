#include "Voronoi.h"

#define JC_VORONOI_IMPLEMENTATION
#include "Utils/jc_voronoi.h"
#define JC_VORONOI_CLIP_IMPLEMENTATION
#include "Utils/jc_voronoi_clip.h"

#include "Utils/Collisions.h"


Voronoi::Voronoi() : Voronoi(std::vector<Vector3>({}))
{

}

Voronoi::Voronoi(int numberRandomPoints, const Vector3& maxBoundarie)
    : Voronoi(numberRandomPoints, Vector3(), maxBoundarie)
{

}

Voronoi::Voronoi(int numberRandomPoints, const Vector3& minBoundarie, const Vector3& maxBoundarie)
    : Voronoi(numberRandomPoints, ShapeCurve({
                                                 Vector3(minBoundarie.xy()),
                                                 Vector3(maxBoundarie.x, minBoundarie.y),
                                                 Vector3(maxBoundarie.xy()),
                                                 Vector3(minBoundarie.x, maxBoundarie.y),
                                             }))
{

}

Voronoi::Voronoi(int numberRandomPoints, ShapeCurve boundingShape)
    : boundingShape(boundingShape.removeDuplicates())
{
    std::tie(this->minBoundarie, this->maxBoundarie) = boundingShape.AABBox();
    this->pointset = boundingShape/*.shrink(1.f)*/.randomPointsInside(numberRandomPoints);
}

Voronoi::Voronoi(std::vector<Vector3> pointset)
    : pointset(pointset)
{
    minBoundarie = Vector3::max();
    maxBoundarie = Vector3::min();
    for (auto& point : pointset) {
        minBoundarie.x = std::min(point.x, minBoundarie.x);
        minBoundarie.y = std::min(point.y, minBoundarie.y);
        minBoundarie.z = std::min(point.z, minBoundarie.z);
        maxBoundarie.x = std::max(point.x, maxBoundarie.x);
        maxBoundarie.y = std::max(point.y, maxBoundarie.y);
        maxBoundarie.z = std::max(point.z, maxBoundarie.z);
    }
    this->boundingShape = ShapeCurve({
                                         Vector3(minBoundarie.xy()),
                                         Vector3(maxBoundarie.x, minBoundarie.y),
                                         Vector3(maxBoundarie.xy()),
                                         Vector3(minBoundarie.x, maxBoundarie.y),
                                     });
}

Voronoi::Voronoi(std::vector<Vector3> pointset, const Vector3& maxBoundarie)
    : Voronoi(pointset, Vector3(0, 0, 0), maxBoundarie)
{

}

Voronoi::Voronoi(std::vector<Vector3> pointset, const Vector3& minBoundarie, const Vector3& maxBoundarie)
    : Voronoi(pointset, ShapeCurve({
                                       Vector3(minBoundarie.xy()),
                                       Vector3(maxBoundarie.x, minBoundarie.y),
                                       Vector3(maxBoundarie.xy()),
                                       Vector3(minBoundarie.x, maxBoundarie.y),
                                   }))
{

}

Voronoi::Voronoi(std::vector<Vector3> pointset, ShapeCurve boundingShape)
    : pointset(pointset), boundingShape(boundingShape.removeDuplicates())
{
    std::tie(this->minBoundarie, this->maxBoundarie) = boundingShape.AABBox();
}

std::vector<ShapeCurve> Voronoi::solve(bool randomizeUntilAllPointsAreSet, int numberOfRelaxations)
{
    this->boundingShape = boundingShape.removeDuplicates();

    if (boundingShape.points.empty()) {
        boundingShape.points = {
            Vector3(minBoundarie.x, minBoundarie.y, minBoundarie.z),
            Vector3(minBoundarie.x, maxBoundarie.y, minBoundarie.z),
            Vector3(maxBoundarie.x, maxBoundarie.y, minBoundarie.z),
            Vector3(maxBoundarie.x, minBoundarie.y, minBoundarie.z)
        };
    }
    if (pointset.size() == 0) {
        return std::vector<ShapeCurve>();
    } else if (pointset.size() == 1) {
        if (!this->boundingShape.points.empty())
            return std::vector<ShapeCurve>{this->boundingShape};
    }
    jcv_point* points = (jcv_point*)malloc( sizeof(jcv_point) * pointset.size());
    jcv_point* boundingPoints = (jcv_point*)malloc( sizeof(jcv_point) * boundingShape.points.size());
    for (size_t i = 0; i < pointset.size(); i++) {
        points[i].x = pointset[i].x;
        points[i].y = pointset[i].y;
        if (i == 0)
            points[i].weight = 1.f;
        else
            points[i].weight = 1.f;
    }
    for (size_t i = 0; i < boundingShape.points.size(); i++) {
        if (i > 0 && boundingShape.points[i] == boundingShape.points[0]) continue;
        boundingPoints[i].x = boundingShape.points[i].x;
        boundingPoints[i].y = boundingShape.points[i].y;
    }

    jcv_clipper polygonclipper;
    jcv_clipper* clipper = 0;
    polygonclipper.test_fn = jcv_clip_polygon_test_point;
    polygonclipper.clip_fn = jcv_clip_polygon_clip_edge;
    polygonclipper.fill_fn = jcv_clip_polygon_fill_gaps;
    jcv_clipping_polygon polygon;
    if (!boundingShape.points.empty()) {

        polygon.num_points = boundingShape.points.size();
        polygon.points = boundingPoints;
    } else {
        std::cout << "We should not be here" << std::endl;
    }
    polygonclipper.ctx = &polygon;

    clipper = &polygonclipper;

    jcv_rect* rect = new jcv_rect;
    rect->min = {minBoundarie.x, minBoundarie.y, 1.f};
    rect->max = {maxBoundarie.x, maxBoundarie.y, 1.f};

    jcv_diagram diagram;
    memset(&diagram, 0, sizeof(jcv_diagram));
    jcv_diagram_generate(pointset.size(), (const jcv_point*)points, rect, clipper, &diagram);

    if (diagram.numsites != (int)pointset.size() && randomizeUntilAllPointsAreSet && numberOfRelaxations > 0) {
        // Let's try new initial points to find at least one configuration with all regions...
        this->pointset = this->boundingShape.randomPointsInside(this->pointset.size());
        jcv_diagram_free(&diagram);
        free(points);
        free(boundingPoints);
        delete rect;
        return this->solve(randomizeUntilAllPointsAreSet, numberOfRelaxations - 1);
    }

    const jcv_site* sites = jcv_diagram_get_sites( &diagram );
    this->areas = std::vector<ShapeCurve>(pointset.size());
    this->neighbors = std::vector<std::vector<int>>(pointset.size());
    for( int i = 0; i < diagram.numsites; ++i )
    {
        const jcv_site* site = &sites[i];
        BSpline areaShape;

        const jcv_graphedge* e = site->edges;
        while( e )
        {
            jcv_point p0 = e->pos[0];
            Vector3 v0 = Vector3(p0.x, p0.y);
            areaShape.points.push_back(v0);
            if (e->neighbor != nullptr)
                neighbors[site->index].push_back(e->neighbor->index);

            e = e->next;
        }
//        if (areaShape.points.size() > 0 && areaShape)
        areas[site->index] = areaShape;
    }

    jcv_diagram_free(&diagram);
    free(points);
    free(boundingPoints);
    delete rect;

    return this->relax(numberOfRelaxations - 1);

}

std::vector<ShapeCurve> Voronoi::relax(int numberOfRelaxations)
{
    for (int relaxation = 0; relaxation < numberOfRelaxations; relaxation++) {
        for (size_t i = 0; i < this->pointset.size(); i++) {
            this->pointset[i] = this->areas[i].centroid();
        }
        this->solve(false, 0);
    }
    return this->areas;
}
