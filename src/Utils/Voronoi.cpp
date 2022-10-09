#include "Voronoi.h"

/*
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/polygon/polygon.hpp>
#include <boost/polygon/voronoi.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/linestring.hpp>
*/

#define JC_VORONOI_IMPLEMENTATION
#include "Utils/jc_voronoi.h"
#define JC_VORONOI_CLIP_IMPLEMENTATION
#include "Utils/jc_voronoi_clip.h"

#include "Utils/Collisions.h"

//using namespace boost::polygon;

Voronoi::Voronoi() : Voronoi(std::vector<Vector3>({}))
{

}

Voronoi::Voronoi(int numberRandomPoints, Vector3 maxBoundarie)
    : Voronoi(numberRandomPoints, Vector3(), maxBoundarie)
{

}

Voronoi::Voronoi(int numberRandomPoints, Vector3 minBoundarie, Vector3 maxBoundarie)
    : Voronoi(numberRandomPoints, ShapeCurve({
                                                 Vector3(minBoundarie.xy()),
                                                 Vector3(minBoundarie.x, maxBoundarie.y),
                                                 Vector3(maxBoundarie.xy()),
                                                 Vector3(maxBoundarie.x, minBoundarie.y)
                                             }))//minBoundarie(minBoundarie), maxBoundarie(maxBoundarie)
{
    /*for (int i = 0; i < numberRandomPoints; i++) {
        pointset.push_back(Vector3(random_gen::generate(minBoundarie.x, maxBoundarie.x),
                                random_gen::generate(minBoundarie.y, maxBoundarie.y)));
    }
    this->boundingShape = ShapeCurve({
                                         Vector3(minBoundarie.xy()),
                                         Vector3(minBoundarie.x, maxBoundarie.y),
                                         Vector3(maxBoundarie.xy()),
                                         Vector3(maxBoundarie.x, minBoundarie.y)
                                     });*/
}

Voronoi::Voronoi(int numberRandomPoints, ShapeCurve boundingShape)
    : boundingShape(boundingShape.removeDuplicates())
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
    this->boundingShape = ShapeCurve({
                                         Vector3(minBoundarie.xy()),
                                         Vector3(minBoundarie.x, maxBoundarie.y),
                                         Vector3(maxBoundarie.xy()),
                                         Vector3(maxBoundarie.x, minBoundarie.y)
                                     });
}

Voronoi::Voronoi(std::vector<Vector3> pointset, Vector3 maxBoundarie)
    : Voronoi(pointset, Vector3(0, 0, 0), maxBoundarie)
{
//    this->maxBoundarie = maxBoundarie;
}

Voronoi::Voronoi(std::vector<Vector3> pointset, Vector3 minBoundarie, Vector3 maxBoundarie)
    : Voronoi(pointset, ShapeCurve({
                                       Vector3(minBoundarie.xy()),
                                       Vector3(minBoundarie.x, maxBoundarie.y),
                                       Vector3(maxBoundarie.xy()),
                                       Vector3(maxBoundarie.x, minBoundarie.y)
                                   }))
{
//    this->minBoundarie = minBoundarie;
//    this->maxBoundarie = maxBoundarie;
}

Voronoi::Voronoi(std::vector<Vector3> pointset, ShapeCurve boundingShape)
    : pointset(pointset), boundingShape(boundingShape.removeDuplicates())
{
    std::tie(this->minBoundarie, this->maxBoundarie) = boundingShape.AABBox();
}

/*
template <>
struct boost::polygon::geometry_concept<Vector3> { typedef point_concept type; };

template <>
struct boost::polygon::point_traits<Vector3> {
  typedef int coordinate_type;

  static inline coordinate_type get(const Vector3& point, boost::polygon::orientation_2d orient) {
    return (orient == HORIZONTAL) ? point.x : point.y;
  }
};
struct Segment {
  Vector3 p0;
  Vector3 p1;
  Segment(Vector3 p0, Vector3 p1) : p0(p0), p1(p1) {}
  Segment(int x1, int y1, int x2, int y2) : p0(x1, y1), p1(x2, y2) {}
};
template <>
struct boost::polygon::geometry_concept<Segment> {
  typedef segment_concept type;
};

template <>
struct boost::polygon::segment_traits<Segment> {
  typedef int coordinate_type;
  typedef Vector3 point_type;

  static inline point_type get(const Segment& segment, direction_1d dir) {
    return dir.to_int() ? segment.p1 : segment.p0;
  }
};

typedef double coordinate_type;
typedef boost::polygon::point_data<coordinate_type> point_type;
typedef boost::polygon::segment_data<coordinate_type> segment_type;
typedef rectangle_data<coordinate_type> rect_type;
typedef voronoi_builder<int> VB;
typedef voronoi_diagram<coordinate_type> VD;
typedef VD::cell_type cell_type;
typedef VD::cell_type::source_index_type source_index_type;
typedef VD::cell_type::source_category_type source_category_type;
typedef VD::edge_type edge_type;
typedef VD::cell_container_type cell_container_type;
typedef VD::cell_container_type vertex_container_type;
typedef VD::edge_container_type edge_container_type;
typedef VD::const_cell_iterator const_cell_iterator;
typedef VD::const_vertex_iterator const_vertex_iterator;
typedef VD::const_edge_iterator const_edge_iterator;

//I'm lazy and use the stl everywhere to avoid writing my own classes
//my toy polygon is a std::vector<Vector3>
typedef std::vector<Vector3> CPolygon;

//we need to specialize our polygon concept mapping in boost polygon
namespace boost { namespace polygon {
  //first register CPolygon as a polygon_concept type
  template <>
  struct geometry_concept<CPolygon>{ typedef polygon_concept type; };

  template <>
  struct polygon_traits<CPolygon> {
    typedef int coordinate_type;
    typedef CPolygon::const_iterator iterator_type;
    typedef Vector3 point_type;

    // Get the begin iterator
    static inline iterator_type begin_points(const CPolygon& t) {
      return t.begin();
    }

    // Get the end iterator
    static inline iterator_type end_points(const CPolygon& t) {
      return t.end();
    }

    // Get the number of sides of the polygon
    static inline std::size_t size(const CPolygon& t) {
      return t.size();
    }

    // Get the winding direction of the polygon
    static inline winding_direction winding(const CPolygon& t) {
      return unknown_winding;
    }
  };

  template <>
  struct polygon_mutable_traits<CPolygon> {
    //expects stl style iterators
    template <typename iT>
    static inline CPolygon& set_points(CPolygon& t,
                                       iT input_begin, iT input_end) {
      t.clear();
      t.insert(t.end(), input_begin, input_end);
      return t;
    }

  };
} }
BOOST_GEOMETRY_REGISTER_POINT_2D(Vector3, double, boost::geometry::cs::cartesian, x, y)
//BOOST_GEOMETRY_REGISTER_POINT_3D(Vector3, double, boost::geometry::cs::cartesian, x, y, z)
BOOST_GEOMETRY_REGISTER_LINESTRING(CPolygon)

point_type retrieve_point(const cell_type& cell, Voronoi& diagram) {
    source_index_type index = cell.source_index();
    return point_data<coordinate_type>(diagram.pointset[index].x, diagram.pointset[index].y);
}

void clip_infinite_edge(const edge_type& edge, std::vector<Vector3>* clipped_edge, Voronoi& diagram) {
    const cell_type& cell1 = *edge.cell();
    const cell_type& cell2 = *edge.twin()->cell();
    point_type origin, direction;

    point_type p1 = retrieve_point(cell1, diagram);
    point_type p2 = retrieve_point(cell2, diagram);
    origin.x((p1.x() + p2.x()) * 0.5);
    origin.y((p1.y() + p2.y()) * 0.5);
    direction.x(p1.y() - p2.y());
    direction.y(p2.x() - p1.x());

    Vector3 vOrigin(origin.x(), origin.y());
    Vector3 vDir = Vector3(direction.x(), direction.y()).normalize();

    Segment s = Segment((edge.vertex0() != NULL ? Vector3(edge.vertex0()->x(), edge.vertex0()->y()) : vOrigin + (vDir * 1000.f)),
                        (edge.vertex1() != NULL ? Vector3(edge.vertex1()->x(), edge.vertex1()->y()) : vOrigin - (vDir * 1000.f)));
    clipped_edge->clear();
    clipped_edge->push_back(s.p0);
    clipped_edge->push_back(s.p1);
}
*/
std::vector<BSpline> Voronoi::solve(bool randomizeUntilAllPointsAreSet, int numberOfRelaxations)
{
    this->boundingShape = boundingShape.removeDuplicates();
    /*
    // Preparing Input Geometries.
    std::vector<Vector3> points = this->pointset;

    // Construction of the Voronoi Diagram.
    boost::polygon::voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), &vd);

    this->areas = std::vector<BSpline>(pointset.size());
    this->neighbors = std::vector<std::vector<int>>(pointset.size());

    unsigned int cell_index = 0;
    for (boost::polygon::voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin(); it != vd.cells().end(); ++it) {
        BSpline area;
        auto e = it->incident_edge();
        if (!e) continue;
        std::cout << it->source_index() << std::endl;
        do {
            std::vector<Vector3> vertices;
            clip_infinite_edge(*e, &vertices, *this);
//            Segment line
//            boost::geometry::intersection()
//            if (e->is_primary()) {
//                auto v = *e->vertex0();
//                area.points.push_back(Vector3(v.x(), v.y()));
            area.points.push_back(vertices.front());
                if (e->twin() && e->twin()->cell()) {
                    neighbors[it->source_index()].push_back(e->twin()->cell()->source_index());
                }
//            }
            e = e->rot_next();
        } while (e != it->incident_edge());

        CPolygon out;
        boost::geometry::intersection(area.points, boundingShape.points, out);
        area.points = out;

        this->areas[it->source_index()] = area;
    }

    return this->relax(numberOfRelaxations - 1);
    */







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
        this->solve(false, 0);
    }
    return this->areas;
}
