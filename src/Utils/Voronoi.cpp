#include "Voronoi.h"

#define JC_VORONOI_IMPLEMENTATION
#include "Utils/jc_voronoi.h"
#define JC_VORONOI_CLIP_IMPLEMENTATION
#include "Utils/jc_voronoi_clip.h"

Voronoi::Voronoi() : Voronoi(std::vector<Vector3>({}))
{

}

Voronoi::Voronoi(int numberRandomPoints, Vector3 minBoundarie, Vector3 maxBoundarie)
{
    for (int i = 0; i < numberRandomPoints; i++) {
        pointset.push_back(Vector3(random_gen::generate(minBoundarie.x, maxBoundarie.x),
                                random_gen::generate(minBoundarie.y, maxBoundarie.y)));
    }
}

Voronoi::Voronoi(std::vector<Vector3> pointset)
    : pointset(pointset)
{
}

std::vector<BSpline> Voronoi::solve()
{
    std::vector<BSpline> areas;

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
    if (!boundingShape.points.empty()) {
        polygonclipper.test_fn = jcv_clip_polygon_test_point;
        polygonclipper.clip_fn = jcv_clip_polygon_clip_edge;
        polygonclipper.fill_fn = jcv_clip_polygon_fill_gaps;

        jcv_clipping_polygon polygon;
        polygon.num_points = boundingShape.points.size();
        polygon.points = boundingPoints;
        polygonclipper.ctx = &polygon;

        clipper = &polygonclipper;
    } else {
        // nada
    }

    jcv_diagram diagram;
    memset(&diagram, 0, sizeof(jcv_diagram));
    jcv_diagram_generate(pointset.size(), (const jcv_point*)points, NULL, clipper, &diagram);

    const jcv_site* sites = jcv_diagram_get_sites( &diagram );
    for( int i = 0; i < diagram.numsites; ++i )
    {
        const jcv_site* site = &sites[i];
        BSpline areaShape;

//        jcv_point s = site->p; // remap(&site->p, &diagram.min, &diagram.max, &dimensions );

        const jcv_graphedge* e = site->edges;
//        std::cout << "Area #" << (i + 1) << "\n";
        while( e )
        {
            jcv_point p0 = e->pos[0]; // remap(&e->pos[0], &diagram.min, &diagram.max, &dimensions );
//            jcv_point p1 = e->pos[1]; // remap(&e->pos[1], &diagram.min, &diagram.max, &dimensions );

            Vector3 v0 = Vector3(p0.x, p0.y);
//            Vector3 v1 = Vector3(p1.x, p1.y);

//            std::cout << v0 << " - " << v1 << "\n";
            areaShape.points.push_back(v0);
            e = e->next;
        }
//        std::cout << std::endl;
        areas.push_back(areaShape);
    }

    jcv_diagram_free(&diagram);
    free(points);
    free(boundingPoints);

    return areas;
}
