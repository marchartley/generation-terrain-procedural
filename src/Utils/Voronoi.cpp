#include "Voronoi.h"

#define JC_VORONOI_IMPLEMENTATION
#include "Utils/jc_voronoi.h"
#define JC_VORONOI_CLIP_IMPLEMENTATION
#include "Utils/jc_voronoi_clip.h"

Voronoi::Voronoi() : Voronoi(std::vector<Vector3>({}))
{

}

Voronoi::Voronoi(int numberRandomPoints, Vector3 maxBoundarie)
    : Voronoi(numberRandomPoints, Vector3(), maxBoundarie)
{

}

Voronoi::Voronoi(int numberRandomPoints, Vector3 minBoundarie, Vector3 maxBoundarie)
    : minBoundarie(minBoundarie), maxBoundarie(maxBoundarie)
{
    for (int i = 0; i < numberRandomPoints; i++) {
        pointset.push_back(Vector3(random_gen::generate(minBoundarie.x, maxBoundarie.x),
                                random_gen::generate(minBoundarie.y, maxBoundarie.y)));
    }
}

Voronoi::Voronoi(int numberRandomPoints, ShapeCurve boundingShape)
    : boundingShape(boundingShape)
{
    std::tie(this->minBoundarie, this->maxBoundarie) = boundingShape.AABBox();
    this->pointset = boundingShape.shrink(1.f).randomPointsInside(numberRandomPoints);
}

Voronoi::Voronoi(std::vector<Vector3> pointset)
    : pointset(pointset)
{
    minBoundarie = Vector3::min();
    maxBoundarie = Vector3::max();
    for (auto& point : pointset) {
        minBoundarie.x = std::min(point.x, minBoundarie.x);
        minBoundarie.y = std::min(point.y, minBoundarie.y);
        minBoundarie.z = std::min(point.z, minBoundarie.z);
        maxBoundarie.x = std::max(point.x, maxBoundarie.x);
        maxBoundarie.y = std::max(point.y, maxBoundarie.y);
        maxBoundarie.z = std::max(point.z, maxBoundarie.z);
    }
}

Voronoi::Voronoi(std::vector<Vector3> pointset, Vector3 maxBoundarie)
    : Voronoi(pointset)
{
    this->maxBoundarie = maxBoundarie;
}

Voronoi::Voronoi(std::vector<Vector3> pointset, Vector3 minBoundarie, Vector3 maxBoundarie)
    : Voronoi(pointset)
{
    this->minBoundarie = minBoundarie;
    this->maxBoundarie = maxBoundarie;
}

std::vector<BSpline> Voronoi::solve(int numberOfRelaxations)
{
    if (boundingShape.points.empty()) {
        boundingShape.points = {
            Vector3(minBoundarie.x, minBoundarie.y, minBoundarie.z),
            Vector3(minBoundarie.x, maxBoundarie.y, minBoundarie.z),
            Vector3(maxBoundarie.x, maxBoundarie.y, minBoundarie.z),
            Vector3(maxBoundarie.x, minBoundarie.y, minBoundarie.z)
        };
    }
    if (pointset.size() == 0) {
        return std::vector<BSpline>();
    } else if (pointset.size() == 1) {
        if (!this->boundingShape.points.empty())
            return std::vector<BSpline>{this->boundingShape};
    }
    jcv_point* points = (jcv_point*)malloc( sizeof(jcv_point) * pointset.size());
    jcv_point* boundingPoints = (jcv_point*)malloc( sizeof(jcv_point) * boundingShape.points.size());
    for (size_t i = 0; i < pointset.size(); i++) {
        points[i].x = pointset[i].x;
        points[i].y = pointset[i].y;
    }
    for (size_t i = 0; i < boundingShape.points.size(); i++) {
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
    rect->min = {minBoundarie.x, minBoundarie.y};
    rect->max = {maxBoundarie.x, maxBoundarie.y};

    jcv_diagram diagram;
    memset(&diagram, 0, sizeof(jcv_diagram));
    jcv_diagram_generate(pointset.size(), (const jcv_point*)points, rect, clipper, &diagram);

    const jcv_site* sites = jcv_diagram_get_sites( &diagram );
    this->areas = std::vector<BSpline>(pointset.size());
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
        areas[site->index] = areaShape;
    }

    jcv_diagram_free(&diagram);
    free(points);
    free(boundingPoints);
    delete rect;

    return this->relax(numberOfRelaxations - 1);
}

std::vector<BSpline> Voronoi::relax(int numberOfRelaxations)
{
    for (int relaxation = 0; relaxation < numberOfRelaxations; relaxation++) {
        for (size_t i = 0; i < this->pointset.size(); i++) {
            this->pointset[i] = this->areas[i].center();
        }
        this->solve();
    }
    return this->areas;
}
