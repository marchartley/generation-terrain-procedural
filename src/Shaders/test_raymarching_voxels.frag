#version 430 core

const float FLT_MIN = 1.175494351e-38;
const float FLT_MAX = 3.402823466e+38;
const float PI_F = 3.14159265358979f;

const int TYPE_SPHERE = 0;
const int TYPE_PLANE = 1;
const int TYPE_SURFACE = 2;
const int TYPE_BOX = 3;
const int TYPE_TORUS = 4;
const int TYPE_RING = 5;
const int TYPE_POINT_LIGHT = 6;
const int TYPE_VOXELS = 7;

const int SHADOW_ENABLED = 1;
const int DBG = 0;
const int DBG_FIRST_VALUE = 1;

const int TOTAL_INTERNAL_REFLECTION = 1;
const int DO_FRESNEL = 1;
const int PLANE_ONESIDE = 1;
const int REFLECT_REDUCE_ITERATION = 1;
const int AMBIENT_OCCLUSION = 0;

const float EPSILON = 0.001;

uniform int NB_SAMPLES = 10;

in vec2 v_texCoord;
in vec4 near_4;    //for computing rays in fragment shader
in vec4 far_4;
in vec3 rayDirection;

uniform sampler3D dataFieldTex;
uniform sampler3D matIndicesTex;
uniform sampler3D matHeightsTex;

struct rt_material {
        vec4 color;
        vec4 absorb;

        float diffuse;
        float reflection;
        float refraction;
        int specular;
        float kd;
        float ks;
};

struct rt_sphere {
        rt_material mat;
        vec4 obj;
        vec4 quat_rotation; // rotate normal
        int textureNum;
        bool hollow;
};

struct rt_plane {
        rt_material mat;
        vec4 pos;
        vec4 normal;
};

struct rt_box {
        rt_material mat;
        vec4 quat_rotation;
        vec4 pos;
        vec4 form;
        int textureNum;
};

struct rt_ring {
        rt_material mat;
        vec4 quat_rotation;
        vec4 pos;
        int textureNum;
        float r1; // square of min radius
        float r2; // square of max radius
};

struct rt_surface {
        rt_material mat;
        vec4 quat_rotation;
        vec4 v_min;
        vec4 v_max;
        vec4 pos;
        float a; // x2
        float b; // y2
        float c; // z2
        float d; // z
        float e; // y
        float f; // const
};

struct rt_torus {
        rt_material mat;
        vec4 quat_rotation;
        vec4 pos;
        vec4 form; // x - radius, y - ring thickness
};

struct rt_light_direct {
        vec4 direction;
        vec4 color;

        float intensity;
};

struct rt_light_point {
        vec4 pos; //pos + radius
        vec4 color;
        float intensity;

        float linear_k;
        float quadratic_k;
};

struct rt_scene {
    vec4 camera_pos;
        mat4 proj_matrix;
        mat4 mv_matrix;
        vec4 bg_color;

        int canvas_width;
        int canvas_height;

        int reflect_depth;
};

struct hit_record {
    rt_material mat;
    vec4 normal;
    float bias_mult;
    float alpha;
};

struct hit_sphere {
    int obj_type;
    int obj_index;
};

uniform int SPHERE_SIZE = 2;
uniform int PLANE_SIZE = 0;
uniform int SURFACE_SIZE = 0;
uniform int BOX_SIZE = 0;
uniform int TORUS_SIZE = 0;
uniform int RING_SIZE = 0;
uniform int LIGHT_DIRECT_SIZE = 0;
uniform int LIGHT_POINT_SIZE = 0;
uniform int AMBIENT_COLOR = 0;
uniform int SHADOW_AMBIENT = 0;
uniform int ITERATIONS = 3;

out vec4 FragColor;


const float maxDist = 1000000.0;

const int MAX_SPHERES_SIZE = 10;
const int MAX_PLANES_SIZE = 10;
const int MAX_SURFACES_SIZE = 10;
const int MAX_BOXES_SIZE = 10;
const int MAX_TORUSES_SIZE = 10;
const int MAX_RINGS_SIZE = 10;
const int MAX_LIGHT_POINTS_SIZE = 10;
const int MAX_DIRECT_LIGHTS_SIZE = 10;

uniform rt_scene scene;
uniform rt_sphere spheres[MAX_SPHERES_SIZE];
uniform rt_plane planes[MAX_PLANES_SIZE];
uniform rt_surface surfaces[MAX_SURFACES_SIZE];
uniform rt_box boxes[MAX_BOXES_SIZE];
uniform rt_torus toruses[MAX_TORUSES_SIZE];
uniform rt_ring rings[MAX_RINGS_SIZE];
uniform rt_light_point light_points[MAX_LIGHT_POINTS_SIZE];
uniform rt_light_direct direct_lights[MAX_DIRECT_LIGHTS_SIZE];


uniform float gridStep = 1.0;

uniform float dx = 0.5;
uniform float dy = 0.5;
uniform float dz = 0.5;

/*
layout( std140, binding = 1 ) uniform scene_buf
{
    rt_scene scene;
};

layout( std140, binding = 2 ) uniform spheres_buf
{
        #if SPHERE_SIZE != 0
    rt_sphere spheres[SPHERE_SIZE];
        #else
        rt_sphere spheres[1];
        #endif
};

layout( std140, binding = 3 ) uniform planes_buf
{
        #if PLANE_SIZE != 0
    rt_plane planes[PLANE_SIZE];
        #else
        rt_plane planes[1];
        #endif
};

layout( std140, binding = 4 ) uniform surfaces_buf
{
        #if SURFACE_SIZE != 0
    rt_surface surfaces[SURFACE_SIZE];
        #else
        rt_surface surfaces[1];
        #endif
};

layout( std140, binding = 5 ) uniform boxes_buf
{
        #if BOX_SIZE != 0
    rt_box boxes[BOX_SIZE];
        #else
        rt_box boxes[1];
        #endif
};

layout( std140, binding = 6 ) uniform toruses_buf
{
        #if TORUS_SIZE != 0
    rt_torus toruses[TORUS_SIZE];
        #else
        rt_torus toruses[1];
        #endif
};

layout( std140, binding = 7 ) uniform rings_buf
{
        #if RING_SIZE != 0
    rt_ring rings[RING_SIZE];
        #else
        rt_ring rings[1];
        #endif
};

layout( std140, binding = 8 ) uniform lights_point_buf
{
        #if LIGHT_POINT_SIZE != 0
    rt_light_point lights_point[LIGHT_POINT_SIZE];
        #else
        rt_light_point lights_point[1];
        #endif
};

layout( std140, binding = 9 ) uniform lights_direct_buf
{
        #if MY_LIGHT_DIRECT_SIZE != 0
    rt_light_direct lights_direct[MY_LIGHT_DIRECT_SIZE];
        #else
        rt_light_direct lights_direct[1];
        #endif
};
*/

void _dbg(float value);
void swap(inout float a, inout float b);
bool isBetween(vec3 value, vec3 min, vec3 max);
vec4 quat_conj(vec4 q);
vec4 quat_inv(vec4 q);
vec4 quat_mult(vec4 q1, vec4 q2);
vec3 rotate(vec4 qr, vec3 v);
vec3 getRayDir();
float dot2( in vec2 v );
float dot2( in vec3 v );
float ndot( in vec2 a, in vec2 b );
float sdPlane( vec3 p, vec3 n, float h );
float rayBoxIntersection(in vec3 orig, in vec3 dir, in vec3 minVert, in vec3 maxVert, in bool returnFirstPositiveValue = false);
bool inBox(in vec3 pos, in vec3 minVert, in vec3 maxVert);
vec3 getWorldCoordinates( in ivec3 _gridCoord );
ivec3 getGridCoordinates( in vec4 _P );
vec3 crossProduct( vec3 a, vec3 b );
void getFirstRayVoxelIntersection( in vec3 origin, in vec3 direction, out ivec3 v0, out vec3 t_n);
vec4 getHittingVoxel(vec3 origin, vec3 direction, sampler3D tex);
int getMaterial(in vec3 pos, in sampler3D matIndicesTex, in sampler3D matHeightsTex, in int invisibleMatIndex, in int maxStack);
vec4 getHittingStackedMaterial(in vec3 origin, in vec3 direction, in sampler3D matIndicesTex, in sampler3D matHeightsTex, in int invisibleMatIndex, out vec3 implicitNormal);
float getKernelVolume(in vec3 pos, in float sigma, in sampler3D matIndicesTex, in sampler3D matHeightsTex, in int invisibleMatIndex, in float maxHeight);
float getMaterialInKernelVolume(in vec3 pos, in float sigma, in sampler3D matIndicesTex, in sampler3D matHeightsTex, in int invisibleMatIndex, in float maxHeight);
float evaluateLayerStackImplicitSurface(in vec3 pos, in float sigma, in int invisibleMatIndex, in float maxHeight);
float sdSphere( vec3 p, float s );
float sdBox( vec3 p, vec3 b );

float sdBoxFrame( vec3 p, vec3 b, float e );
float sdEllipsoid( in vec3 p, in vec3 r );
float sdTorus( vec3 p, vec2 t );
float sdCappedTorus(in vec3 p, in vec2 sc, in float ra, in float rb);
float sdHexPrism( vec3 p, vec2 h );
float sdOctogonPrism( in vec3 p, in float r, float h );
float sdCapsule( vec3 p, vec3 a, vec3 b, float r );
float sdRoundCone( in vec3 p, in float r1, float r2, float h );
float sdRoundCone(vec3 p, vec3 a, vec3 b, float r1, float r2);
float sdTriPrism( vec3 p, vec2 h );
float sdCylinder( vec3 p, vec2 h );
float sdCylinder(vec3 p, vec3 a, vec3 b, float r);
float sdCone( in vec3 p, in vec2 c, float h );
float sdCappedCone( in vec3 p, in float h, in float r1, in float r2 );
float sdCappedCone(vec3 p, vec3 a, vec3 b, float ra, float rb);
float sdSolidAngle(vec3 pos, vec2 c, float ra);
float sdOctahedron(vec3 p, float s);
float sdPyramid( in vec3 p, in float h );
float sdRhombus(vec3 p, float la, float lb, float h, float ra);
float sdHorseshoe( in vec3 p, in vec2 c, in float r, in float le, vec2 w );
float sdU( in vec3 p, in float r, in float le, vec2 w );
vec2 iBox( in vec3 ro, in vec3 rd, in vec3 rad );

vec2 opU( vec2 d1, vec2 d2 );
float opUnion( float d1, float d2 );
float opSubtraction( float d1, float d2 );
float opIntersection( float d1, float d2 );
float opSmoothUnion( float d1, float d2, float k );
float opSmoothSubtraction( float d1, float d2, float k );
float opSmoothIntersection( float d1, float d2, float k );

vec2 map( in vec3 pos, in vec3 rd, out hit_sphere hit);
vec3 raycast( in vec3 ro, in vec3 rd , out vec3 normal );
float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax );
vec3 calcNormal( in vec3 pos );
float calcAO( in vec3 pos, in vec3 nor );
vec3 render( in vec3 ro, in vec3 rd, in vec3 rdx, in vec3 rdy );

float checkersGradBox( in vec2 p, in vec2 dpdx, in vec2 dpdy );


void _dbg(float value)
{
    FragColor = vec4(value,value,value,1);
}

void _dbg(vec3 value)
{
    FragColor = vec4(value, 1.0);
}

void swap(inout float a, inout float b)
{
        float tmp = a;
        a = b;
        b = tmp;
}

bool isBetween(vec3 value, vec3 min, vec3 max)
{
        return greaterThan(value, min) == bvec3(true) && lessThan(value, max) == bvec3(true);
}

vec4 quat_conj(vec4 q)
{
        return vec4(-q.x, -q.y, -q.z, q.w);
}

vec4 quat_inv(vec4 q)
{
        return quat_conj(q) * (1 / dot(q, q));
}

vec4 quat_mult(vec4 q1, vec4 q2)
{
        vec4 qr;
        qr.x = (q1.w * q2.x) + (q1.x * q2.w) + (q1.y * q2.z) - (q1.z * q2.y);
        qr.y = (q1.w * q2.y) - (q1.x * q2.z) + (q1.y * q2.w) + (q1.z * q2.x);
        qr.z = (q1.w * q2.z) + (q1.x * q2.y) - (q1.y * q2.x) + (q1.z * q2.w);
        qr.w = (q1.w * q2.w) - (q1.x * q2.x) - (q1.y * q2.y) - (q1.z * q2.z);
        return qr;
}

vec3 rotate(vec4 qr, vec3 v)
{
        vec4 qr_conj = quat_conj(qr);
        vec4 q_pos = vec4(v.xyz, 0);
        vec4 q_tmp = quat_mult(qr, q_pos);
        return quat_mult(q_tmp, qr_conj).xyz;
}

vec3 getRayDir()
{
    return rayDirection;
}

float dot2( in vec2 v ) { return dot(v,v); }
float dot2( in vec3 v ) { return dot(v,v); }
float ndot( in vec2 a, in vec2 b ) { return a.x*b.x - a.y*b.y; }

float sdPlane( vec3 p, vec3 n, float h )
{
  // n must be normalized
  return dot(p, n) + h;
}


float rayBoxIntersection(in vec3 orig, in vec3 dir, in vec3 minVert, in vec3 maxVert, in bool returnFirstPositiveValue = false) {
    float tmin = (minVert.x - orig.x) / dir.x;
    float tmax = (maxVert.x - orig.x) / dir.x;

    if (tmin > tmax) swap(tmin, tmax);

    float tymin = (minVert.y - orig.y) / dir.y;
    float tymax = (maxVert.y - orig.y) / dir.y;

    if (tymin > tymax) swap(tymin, tymax);

    if ((tmin > tymax) || (tymin > tmax))
        return -1;

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    float tzmin = (minVert.z - orig.z) / dir.z;
    float tzmax = (maxVert.z - orig.z) / dir.z;

    if (tzmin > tzmax) swap(tzmin, tzmax);

    if ((tmin > tzmax) || (tzmin > tmax))
        return -1;

    if (tzmin > tmin)
        tmin = tzmin;

    if (tzmax < tmax)
        tmax = tzmax;

    return returnFirstPositiveValue ? (tmin > 0.0 ? tmin : tmax) : min(tmin, tmax);
}

bool inBox(in vec3 pos, in vec3 minVert, in vec3 maxVert) {
    return isBetween(pos, minVert, maxVert);
}

vec3 getWorldCoordinates( in ivec3 _gridCoord ){
    return vec3( (_gridCoord.x+0.5)*dx,
                 (_gridCoord.y+0.5)*dy,
                 (_gridCoord.z+0.5)*dz );
}

ivec3 getGridCoordinates( in vec4 _P ){
    return ivec3( int( _P.x/dx ) ,
                  int( _P.y/dy ) ,
                  int( _P.z/dz ) );
}

vec3 crossProduct( vec3 a, vec3 b ){

    return vec3( a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x );
}

void getFirstRayVoxelIntersection( in vec3 origin, in vec3 direction, out ivec3 v0, out vec3 t_n){
    v0 = getGridCoordinates(vec4(origin.xyz, 1.));

    float xi = v0.x*dx;
    if( direction.x > 0  ) xi = xi + dx;
    float yi = v0.y*dy;
    if( direction.y > 0  ) yi = yi + dy;
    float zi = v0.z*dz;
    if( direction.z > 0  ) zi = zi + dz;
    t_n = vec3 ( ((xi - origin.x)/direction.x), ((yi - origin.y)/direction.y), ((zi - origin.z)/direction.z) );

    if( abs( direction.x ) < 0.00001 ) t_n.x = 100000000;
    if( abs( direction.y ) < 0.00001 ) t_n.y = 100000000;
    if( abs( direction.z ) < 0.00001 ) t_n.z = 100000000;

}

vec4 getHittingVoxel(vec3 origin, vec3 direction, sampler3D tex) {
    vec3 V = direction;

    bool in_tet = true;
    bool hit = false;


    vec3 t_next;
    ivec3 next_voxel;
    ivec3 origin_voxel;

    /**************initialization******************/

    vec3 texSize = textureSize(tex, 0);
    vec3 P = origin + V * (rayBoxIntersection(origin, V, vec3(0.0), texSize, true) - 1.0);
    vec3 Current_P = P.xyz;

    //Find the first intersection of the ray with the grid
    getFirstRayVoxelIntersection(Current_P, V, origin_voxel, t_next );

    vec3 dt = vec3( abs(dx/V.x), abs(dy/V.y), abs(dz/V.z) );

    ivec3 grid_step = ivec3 (-1, -1, -1);


    if( V.x > 0 ) grid_step.x = 1;
    if( V.y > 0 ) grid_step.y = 1;
    if( V.z > 0 ) grid_step.z = 1;

    float t_min = 0;

    /***********************************************/

    vec3 Current_text3DCoord;

    vec4 Pos;

    vec4 color = vec4 (0.,0.,0.,1.);

    float v_step = dx;
    if( v_step > dy ) v_step = dy;
    if( v_step > dz ) v_step = dz;

    bool changed = false;
    next_voxel = origin_voxel;


    vec3 normals[6];
    normals[0] = vec3( 1., 0., 0.); normals[1] = vec3( -1.,  0.,  0.);
    normals[2] = vec3( 0., 1., 0.); normals[3] = vec3(  0., -1.,  0.);
    normals[4] = vec3( 0., 0., 1.); normals[5] = vec3(  0.,  0., -1.);

    vec3 n;

    int normal_index = -1;

    float val =0;
    int fragmentIteration = 0;
    while( in_tet && !hit/* && fragmentIteration < 100 */){
        fragmentIteration++;
        if( t_next.x < t_next.y && t_next.x < t_next.z ){
            Current_P = P.xyz + t_next.x*V;
            t_min = t_next.x;
            t_next.x = t_next.x + dt.x;
            next_voxel.x = next_voxel.x + grid_step.x;
            if( V.x > 0 ) {
                n = normals[1];
                normal_index = 1;
            } else {
                n = normals[0];
                normal_index = 0;
            }
        } else if( t_next.y < t_next.x && t_next.y < t_next.z ){
            Current_P = P.xyz + t_next.y*V;
            t_min = t_next.y;
            t_next.y = t_next.y + dt.y;
            next_voxel.y = next_voxel.y + grid_step.y;
            if( V.y > 0 ) {
                n = normals[3];
                normal_index = 3;
            } else {
                n = normals[2];
                normal_index = 2;
            }
        } else{
            Current_P = P.xyz + t_next.z*V;
            t_min = t_next.z;
            t_next.z = t_next.z + dt.z;
            next_voxel.z = next_voxel.z + grid_step.z;
            if( V.z > 0 ) {
                n = normals[5];
                normal_index = 5;
            } else {
                n = normals[4];
                normal_index = 4;
            }
        }

        Current_P = Current_P + 0.0001*v_step*V;
        vec3 voxel_center_P = getWorldCoordinates( next_voxel );

        if (texelFetch(tex, ivec3(voxel_center_P), 0).a > 0.5) {
            hit = true;
        }

        if (!inBox(Current_P, vec3(0.0) - vec3(5.0), texSize + vec3(5.0))) {
            in_tet = false;
        }
    }

    if (!in_tet || !hit) {
        // = out of box or many iterations without a hit
        return vec4(maxDist, maxDist, maxDist, -1);
    }
    return vec4(Current_P, normal_index);
}

int getMaterial(in vec3 pos, in sampler3D matIndicesTex, in sampler3D matHeightsTex, in int invisibleMatIndex, in int maxStack) {
    float z = pos.z;
    ivec2 xy = ivec2(pos.xy);
    float currentHeight = 0;
    int currentStackedIndex = 0;
    while (currentHeight < z && currentStackedIndex < maxStack) {
        currentHeight += texelFetch(matHeightsTex, ivec3(xy, currentStackedIndex), 0).a;
        currentStackedIndex ++;
    }
    // If a layer was found, return it. Otherwise, return the invisible material
    return (currentStackedIndex < maxStack ? int(texelFetch(matIndicesTex, ivec3(xy, currentStackedIndex - 1), 0).a) : invisibleMatIndex);
}

vec4 getHittingStackedMaterial(in vec3 origin, in vec3 direction, in sampler3D matIndicesTex, in sampler3D matHeightsTex, in int invisibleMatIndex, out vec3 implicitNormal) {
    vec3 V = direction;
    int maxStack = textureSize(matIndicesTex, 0).z;

    bool in_tet = true;
    bool hit = false;


    vec3 t_next;
    vec3 next_voxel;
    vec3 origin_voxel;

    /**************initialization******************/

    vec3 texSize = vec3(textureSize(matIndicesTex, 0).xy, 150); // TODO : Compute the max height of the terrain
    vec3 P = origin + V * (max(0.0, rayBoxIntersection(origin, V, vec3(0.0), texSize, false)) - 1.0);
    vec3 Current_P = P.xyz;

    //Find the first intersection of the ray with the grid
//    getFirstRayVoxelIntersection(Current_P, V, origin_voxel, t_next );
    origin_voxel = Current_P; t_next = V;

    vec3 dt = vec3( abs(dx/V.x), abs(dy/V.y), abs(dz/V.z) );

    ivec3 grid_step = ivec3 (-1, -1, -1);


    if( V.x > 0 ) grid_step.x = 1;
    if( V.y > 0 ) grid_step.y = 1;
    if( V.z > 0 ) grid_step.z = 1;

    float t_min = 0;

    /***********************************************/

    vec3 Current_text3DCoord;

    vec4 Pos;

    vec4 color = vec4 (0.,0.,0.,1.);

    float v_step = dx;
    if( v_step > dy ) v_step = dy;
    if( v_step > dz ) v_step = dz;

    bool changed = false;
    next_voxel = origin_voxel;


    vec3 normals[6];
    normals[0] = vec3( 1., 0., 0.); normals[1] = vec3( -1.,  0.,  0.);
    normals[2] = vec3( 0., 1., 0.); normals[3] = vec3(  0., -1.,  0.);
    normals[4] = vec3( 0., 0., 1.); normals[5] = vec3(  0.,  0., -1.);

    vec3 n;

    int normal_index = -1;

    float val =0;
    int fragmentIteration = 0;
    while( in_tet && !hit/* && fragmentIteration < 100 */){
        fragmentIteration++;
        if( t_next.x < t_next.y && t_next.x < t_next.z ){
            Current_P = P.xyz + t_next.x*V;
            t_min = t_next.x;
            t_next.x = t_next.x + dt.x;
            next_voxel.x = next_voxel.x + grid_step.x;
            if( V.x > 0 ) {
                n = normals[1];
                normal_index = 1;
            } else {
                n = normals[0];
                normal_index = 0;
            }
        } else if( t_next.y < t_next.x && t_next.y < t_next.z ){
            Current_P = P.xyz + t_next.y*V;
            t_min = t_next.y;
            t_next.y = t_next.y + dt.y;
            next_voxel.y = next_voxel.y + grid_step.y;
            if( V.y > 0 ) {
                n = normals[3];
                normal_index = 3;
            } else {
                n = normals[2];
                normal_index = 2;
            }
        } else{
            Current_P = P.xyz + t_next.z*V;
            t_min = t_next.z;
            t_next.z = t_next.z + dt.z;
            next_voxel.z = next_voxel.z + grid_step.z;
            if( V.z > 0 ) {
                n = normals[5];
                normal_index = 5;
            } else {
                n = normals[4];
                normal_index = 4;
            }
        }

//        Current_P = Current_P + 0.0001*v_step*V;
        Current_P += V;
//        vec3 voxel_center_P = getWorldCoordinates( next_voxel );

        float implicit = evaluateLayerStackImplicitSurface(Current_P, 1.0, 0, texSize.z);
//        if (getMaterial(Current_P, matIndicesTex, matHeightsTex, invisibleMatIndex, maxStack) != invisibleMatIndex) {
        if (abs(implicit) <= 0.5) {
            hit = true;
        }

        if (!inBox(Current_P, vec3(0.0) - vec3(5.0), texSize + vec3(5.0))) {
            in_tet = false;
        }
    }

    if (!in_tet || !hit) {
        // = out of box or many iterations without a hit
        return vec4(maxDist, maxDist, maxDist, -1);
    }
    implicitNormal = vec3(
                evaluateLayerStackImplicitSurface(Current_P + vec3(0.01, 0, 0), 1.0, 0, texSize.z) - evaluateLayerStackImplicitSurface(Current_P - vec3(0.01, 0, 0), 1.0, 0, texSize.z),
                evaluateLayerStackImplicitSurface(Current_P + vec3(0, 0.01, 0), 1.0, 0, texSize.z) - evaluateLayerStackImplicitSurface(Current_P - vec3(0, 0.01, 0), 1.0, 0, texSize.z),
                evaluateLayerStackImplicitSurface(Current_P + vec3(0, 0, 0.01), 1.0, 0, texSize.z) - evaluateLayerStackImplicitSurface(Current_P - vec3(0, 0, 0.01), 1.0, 0, texSize.z)
                );
    return vec4(Current_P, normal_index);
}

float getKernelVolume(in vec3 pos, in float sigma, in sampler3D matIndicesTex, in sampler3D matHeightsTex, in int invisibleMatIndex, in float maxHeight)
{
    // Most of the parameters are useless, but it's for consistency with next function
//    vec3 minVert = min(pos - vec3(sigma), vec3(0.0));
//    vec3 maxVert = max(pos + vec3(sigma), vec3(textureSize(matIndicesTex, 0).xy, maxHeight));

    return 8 * sigma * sigma * sigma;
}

float getMaterialInKernelVolume(in vec3 pos, in float sigma, in sampler3D matIndicesTex, in sampler3D matHeightsTex, in int invisibleMatIndex, in float maxHeight)
{
    // Most of the parameters are useless, but it's for consistency with next function
    vec3 texSize = textureSize(matIndicesTex, 0).xyz;
    float totalVolume = 0.0;
    float minX = pos.x - sigma; // max(pos.x - sigma, 0.0);
    float maxX = pos.x + sigma; // min(pos.x + sigma, texSize.x);
    float minY = pos.y - sigma; // max(pos.y - sigma, 0.0);
    float maxY = pos.y + sigma; // min(pos.y + sigma, texSize.y);
    float minZ = pos.z - sigma; // max(pos.z - sigma, 0.0);
    float maxZ = pos.z + sigma;

    float count = 0.0;
    for (int x = int(floor(minX)); x < int(ceil(maxX)); x++) {
        for (int y = int(floor(minY)); y < int(ceil(maxY)); y++) {
            count += 1.0;

            float widthX = 1.0; // min(x+1, maxX) - max(x, minX);
            float widthY = 1.0; // min(y+1, maxY) - max(y, minY);

            float currentHeight = 0.0;
            for (int i = 0; i < texSize.z/* && currentHeight < maxZ*/; i++) {
                float blockHeight = texelFetch(matHeightsTex, ivec3(x, y, i), 0).a;
                int blockMaterial = int(texelFetch(matIndicesTex, ivec3(x, y, i), 0).a);

                float start = clamp(currentHeight, minZ, maxZ);
                float end   = clamp(currentHeight + blockHeight, minZ, maxZ);
                totalVolume += (end - start) * (blockMaterial == invisibleMatIndex ? 0.0 : 1.0) * (widthX * widthY);
                currentHeight += blockHeight;
            }
        }
    }
//    return count * (2 * sigma);
    return totalVolume;
}

float evaluateLayerStackImplicitSurface(in vec3 pos, in float sigma, in int invisibleMatIndex, in float maxHeight)
{
    float V_all = getKernelVolume(pos, sigma, matIndicesTex, matHeightsTex, invisibleMatIndex, maxHeight);
    float V_material = getMaterialInKernelVolume(pos, sigma, matIndicesTex, matHeightsTex, invisibleMatIndex, maxHeight);

    return 2 * (V_material / V_all) - 1.0; // Gives value between -1 (empty) and 1 (full)
}


float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

float sdBox( vec3 p, vec3 b )
{
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float sdBoxFrame( vec3 p, vec3 b, float e )
{
       p = abs(p  )-b;
  vec3 q = abs(p+e)-e;

  return min(min(
      length(max(vec3(p.x,q.y,q.z),0.0))+min(max(p.x,max(q.y,q.z)),0.0),
      length(max(vec3(q.x,p.y,q.z),0.0))+min(max(q.x,max(p.y,q.z)),0.0)),
      length(max(vec3(q.x,q.y,p.z),0.0))+min(max(q.x,max(q.y,p.z)),0.0));
}

float sdEllipsoid( in vec3 p, in vec3 r ) // approximated
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float sdTorus( vec3 p, vec2 t )
{
    return length( vec2(length(p.xz)-t.x,p.y) )-t.y;
}

float sdCappedTorus(in vec3 p, in vec2 sc, in float ra, in float rb)
{
    p.x = abs(p.x);
    float k = (sc.y*p.x>sc.x*p.y) ? dot(p.xy,sc) : length(p.xy);
    return sqrt( dot(p,p) + ra*ra - 2.0*ra*k ) - rb;
}

float sdHexPrism( vec3 p, vec2 h )
{
    vec3 q = abs(p);

    const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
    p = abs(p);
    p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
    vec2 d = vec2(
       length(p.xy - vec2(clamp(p.x, -k.z*h.x, k.z*h.x), h.x))*sign(p.y - h.x),
       p.z-h.y );
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdOctogonPrism( in vec3 p, in float r, float h )
{
  const vec3 k = vec3(-0.9238795325,   // sqrt(2+sqrt(2))/2
                       0.3826834323,   // sqrt(2-sqrt(2))/2
                       0.4142135623 ); // sqrt(2)-1
  // reflections
  p = abs(p);
  p.xy -= 2.0*min(dot(vec2( k.x,k.y),p.xy),0.0)*vec2( k.x,k.y);
  p.xy -= 2.0*min(dot(vec2(-k.x,k.y),p.xy),0.0)*vec2(-k.x,k.y);
  // polygon side
  p.xy -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
  vec2 d = vec2( length(p.xy)*sign(p.y), p.z-h );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
        vec3 pa = p-a, ba = b-a;
        float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
        return length( pa - ba*h ) - r;
}

float sdRoundCone( in vec3 p, in float r1, float r2, float h )
{
    vec2 q = vec2( length(p.xz), p.y );

    float b = (r1-r2)/h;
    float a = sqrt(1.0-b*b);
    float k = dot(q,vec2(-b,a));

    if( k < 0.0 ) return length(q) - r1;
    if( k > a*h ) return length(q-vec2(0.0,h)) - r2;

    return dot(q, vec2(a,b) ) - r1;
}

float sdRoundCone(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    // sampling independent computations (only depend on shape)
    vec3  ba = b - a;
    float l2 = dot(ba,ba);
    float rr = r1 - r2;
    float a2 = l2 - rr*rr;
    float il2 = 1.0/l2;

    // sampling dependant computations
    vec3 pa = p - a;
    float y = dot(pa,ba);
    float z = y - l2;
    float x2 = dot2( pa*l2 - ba*y );
    float y2 = y*y*l2;
    float z2 = z*z*l2;

    // single square root!
    float k = sign(rr)*rr*rr*x2;
    if( sign(z)*a2*z2 > k ) return  sqrt(x2 + z2)        *il2 - r2;
    if( sign(y)*a2*y2 < k ) return  sqrt(x2 + y2)        *il2 - r1;
                            return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}

float sdTriPrism( vec3 p, vec2 h )
{
    const float k = sqrt(3.0);
    h.x *= 0.5*k;
    p.xy /= h.x;
    p.x = abs(p.x) - 1.0;
    p.y = p.y + 1.0/k;
    if( p.x+k*p.y>0.0 ) p.xy=vec2(p.x-k*p.y,-k*p.x-p.y)/2.0;
    p.x -= clamp( p.x, -2.0, 0.0 );
    float d1 = length(p.xy)*sign(-p.y)*h.x;
    float d2 = abs(p.z)-h.y;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

// vertical
float sdCylinder( vec3 p, vec2 h )
{
    vec2 d = abs(vec2(length(p.xz),p.y)) - h;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

// arbitrary orientation
float sdCylinder(vec3 p, vec3 a, vec3 b, float r)
{
    vec3 pa = p - a;
    vec3 ba = b - a;
    float baba = dot(ba,ba);
    float paba = dot(pa,ba);

    float x = length(pa*baba-ba*paba) - r*baba;
    float y = abs(paba-baba*0.5)-baba*0.5;
    float x2 = x*x;
    float y2 = y*y*baba;
    float d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    return sign(d)*sqrt(abs(d))/baba;
}

// vertical
float sdCone( in vec3 p, in vec2 c, float h )
{
    vec2 q = h*vec2(c.x,-c.y)/c.y;
    vec2 w = vec2( length(p.xz), p.y );

        vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
    vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
    float k = sign( q.y );
    float d = min(dot( a, a ),dot(b, b));
    float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
        return sqrt(d)*sign(s);
}

float sdCappedCone( in vec3 p, in float h, in float r1, in float r2 )
{
    vec2 q = vec2( length(p.xz), p.y );

    vec2 k1 = vec2(r2,h);
    vec2 k2 = vec2(r2-r1,2.0*h);
    vec2 ca = vec2(q.x-min(q.x,(q.y < 0.0)?r1:r2), abs(q.y)-h);
    vec2 cb = q - k1 + k2*clamp( dot(k1-q,k2)/dot2(k2), 0.0, 1.0 );
    float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
    return s*sqrt( min(dot2(ca),dot2(cb)) );
}

float sdCappedCone(vec3 p, vec3 a, vec3 b, float ra, float rb)
{
    float rba  = rb-ra;
    float baba = dot(b-a,b-a);
    float papa = dot(p-a,p-a);
    float paba = dot(p-a,b-a)/baba;

    float x = sqrt( papa - paba*paba*baba );

    float cax = max(0.0,x-((paba<0.5)?ra:rb));
    float cay = abs(paba-0.5)-0.5;

    float k = rba*rba + baba;
    float f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );

    float cbx = x-ra - f*rba;
    float cby = paba - f;

    float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;

    return s*sqrt( min(cax*cax + cay*cay*baba,
                       cbx*cbx + cby*cby*baba) );
}

// c is the sin/cos of the desired cone angle
float sdSolidAngle(vec3 pos, vec2 c, float ra)
{
    vec2 p = vec2( length(pos.xz), pos.y );
    float l = length(p) - ra;
        float m = length(p - c*clamp(dot(p,c),0.0,ra) );
    return max(l,m*sign(c.y*p.x-c.x*p.y));
}

float sdOctahedron(vec3 p, float s)
{
    p = abs(p);
    float m = p.x + p.y + p.z - s;

    // exact distance
    #if 0
    vec3 o = min(3.0*p - m, 0.0);
    o = max(6.0*p - m*2.0 - o*3.0 + (o.x+o.y+o.z), 0.0);
    return length(p - s*o/(o.x+o.y+o.z));
    #endif

    // exact distance
    #if 1
        vec3 q;
         if( 3.0*p.x < m ) q = p.xyz;
    else if( 3.0*p.y < m ) q = p.yzx;
    else if( 3.0*p.z < m ) q = p.zxy;
    else return m*0.57735027;
    float k = clamp(0.5*(q.z-q.y+s),0.0,s);
    return length(vec3(q.x,q.y-s+k,q.z-k));
    #endif

    // bound, not exact
    #if 0
        return m*0.57735027;
    #endif
}

float sdPyramid( in vec3 p, in float h )
{
    float m2 = h*h + 0.25;

    // symmetry
    p.xz = abs(p.xz);
    p.xz = (p.z>p.x) ? p.zx : p.xz;
    p.xz -= 0.5;

    // project into face plane (2D)
    vec3 q = vec3( p.z, h*p.y - 0.5*p.x, h*p.x + 0.5*p.y);

    float s = max(-q.x,0.0);
    float t = clamp( (q.y-0.5*p.z)/(m2+0.25), 0.0, 1.0 );

    float a = m2*(q.x+s)*(q.x+s) + q.y*q.y;
        float b = m2*(q.x+0.5*t)*(q.x+0.5*t) + (q.y-m2*t)*(q.y-m2*t);

    float d2 = min(q.y,-q.x*m2-q.y*0.5) > 0.0 ? 0.0 : min(a,b);

    // recover 3D and scale, and add sign
    return sqrt( (d2+q.z*q.z)/m2 ) * sign(max(q.z,-p.y));;
}

// la,lb=semi axis, h=height, ra=corner
float sdRhombus(vec3 p, float la, float lb, float h, float ra)
{
    p = abs(p);
    vec2 b = vec2(la,lb);
    float f = clamp( (ndot(b,b-2.0*p.xz))/dot(b,b), -1.0, 1.0 );
        vec2 q = vec2(length(p.xz-0.5*b*vec2(1.0-f,1.0+f))*sign(p.x*b.y+p.z*b.x-b.x*b.y)-ra, p.y-h);
    return min(max(q.x,q.y),0.0) + length(max(q,0.0));
}

float sdHorseshoe( in vec3 p, in vec2 c, in float r, in float le, vec2 w )
{
    p.x = abs(p.x);
    float l = length(p.xy);
    p.xy = mat2(-c.x, c.y,
              c.y, c.x)*p.xy;
    p.xy = vec2((p.y>0.0 || p.x>0.0)?p.x:l*sign(-c.x),
                (p.x>0.0)?p.y:l );
    p.xy = vec2(p.x,abs(p.y-r))-vec2(le,0.0);

    vec2 q = vec2(length(max(p.xy,0.0)) + min(0.0,max(p.x,p.y)),p.z);
    vec2 d = abs(q) - w;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdU( in vec3 p, in float r, in float le, vec2 w )
{
    p.x = (p.y>0.0) ? abs(p.x) : length(p.xy);
    p.x = abs(p.x-r);
    p.y = p.y - le;
    float k = max(p.x,p.y);
    vec2 q = vec2( (k<0.0) ? -k : length(max(p.xy,0.0)), abs(p.z) ) - w;
    return length(max(q,0.0)) + min(max(q.x,q.y),0.0);
}

//------------------------------------------------------------------

vec2 opU( vec2 d1, vec2 d2 )
{
        return (d1.x<d2.x) ? d1 : d2;
}

float opUnion( float d1, float d2 ) { return min(d1,d2); }

float opSubtraction( float d1, float d2 ) { return max(-d1,d2); }

float opIntersection( float d1, float d2 ) { return max(d1,d2); }

float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); }

float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h); }

float opSmoothIntersection( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) + k*h*(1.0-h); }
//------------------------------------------------------------------

vec2 map( in vec3 pos, in vec3 rd, out hit_sphere hit)
{
    float minDist = maxDist;
    hit.obj_index = -1;
    vec2 res = vec2( minDist, 0.0 );
    float dist;

    for (int i = 0; i < SPHERE_SIZE; i++) {
        dist = sdSphere(pos - spheres[i].obj.xyz, spheres[i].obj.w );
        minDist = opUnion(minDist, dist);

        if (dist > 0 && dist < EPSILON) {
            hit.obj_type = TYPE_SPHERE;
            hit.obj_index = i;
        }
    }
    for (int i = 0; i < PLANE_SIZE; i++) {
        dist = sdPlane(pos - planes[i].pos.xyz, planes[i].normal.xyz, 1.0); // ARGUMENT H IS WRONG!!
        minDist = opUnion(minDist, dist);

        if (dist > 0 && dist < EPSILON) {
            hit.obj_type = TYPE_PLANE;
            hit.obj_index = i;
        }
    }
    for (int i = 0; i < SURFACE_SIZE; i++) {
//        dist = sd..(); // NEED TO DEFINE A SURFACE A LITTLE BETTER
//        minDist = opUnion(minDist, dist);

//        if (dist > 0 && dist < EPSILON) {
//            hit.obj_type = TYPE_SURFACE;
//            hit.obj_index = i;
//        }
    }
    for (int i = 0; i < BOX_SIZE; i++) {
        dist = sdBox(pos - boxes[i].pos.xyz, boxes[i].form.xyz); // UPDATE ROTATION
        minDist = opUnion(minDist, dist);

        if (dist > 0 && dist < EPSILON) {
            hit.obj_type = TYPE_BOX;
            hit.obj_index = i;
        }
    }
    for (int i = 0; i < TORUS_SIZE; i++) {
        dist = sdTorus(pos - toruses[i].pos.xyz, toruses[i].form.xy);
        minDist = opUnion(minDist, dist);

        if (dist > 0 && dist < EPSILON) {
            hit.obj_type = TYPE_TORUS;
            hit.obj_index = i;
        }
    }
    for (int i = 0; i < RING_SIZE; i++) {
//        dist = sdRing(); // FIND FUNCTION
//        minDist = opUnion(minDist, dist);

//        if (dist > 0 && dist < EPSILON) {
//            hit.obj_type = TYPE_RING;
//            hit.obj_index = i;
//        }
    }

    return vec2(minDist, 5.0);

//    // Voxels
//    vec3 texSize = textureSize(dataFieldTex, 0);

////    dist = sdVoxels(pos, rd, dataFieldTex);
//    vec4 hittingVoxel = getHittingVoxel(pos, rd, dataFieldTex);
//    dist = length(hittingVoxel.xyz - pos);
//    minDist = opUnion(minDist, dist);

//    if (dist > 0 && dist < EPSILON) {
//        hit.obj_type = TYPE_VOXELS;
//        hit.obj_index = 0;
//    }
//    res = opU(res, vec2( sdBox(pos - texSize * .5, texSize * .5), -10.0));
//    for (int x = -1; x <= 1; x++)
//        for (int y = -1; y <= 1; y++)
//            for (int z = -1; z <= 1; z++)
//                if (inBox(pos + vec3(x, y, z), vec3(0.0), texSize) && texelFetch(dataFieldTex, ivec3(pos + vec3(x, y, z)), 0).a > .5)
//                    res = opU( res, vec2(sdBox((pos - ivec3(pos + vec3(x, y, z))), vec3(.5)), 5.0));


//    float resolution = 0.6;
//    for (float x = -1.0; x < 1.0; x += resolution)
//        for (float y = -1.0; y < 1.0; y += resolution)
//            for (float z = -1.0; z < 1.0; z += resolution) {
//                vec3 checkPos = pos + vec3(x, y, z);
//                if (texture(dataFieldTex, checkPos / texSize).a > 0.5) {
//                    res = opU(res, vec2( sdSphere(vec3(x, y, z), 1.0), 10.0 ));
//                }
//            }

    return res;
/*
    // bounding box
    if( sdBox( pos-vec3(-2.0,0.3,0.25),vec3(0.3,0.3,1.0) )<res.x )
    {
      res = opU( res, vec2( sdSphere(    pos-vec3(-2.0,0.25, 0.0), 0.25 ), 26.9 ) );
          res = opU( res, vec2( sdRhombus(  (pos-vec3(-2.0,0.25, 1.0)).xzy, 0.15, 0.25, 0.04, 0.08 ),17.0 ) );
    }

    // bounding box
    if( sdBox( pos-vec3(0.0,0.3,-1.0),vec3(0.35,0.3,2.5) )<res.x )
    {
        res = opU( res, vec2( sdCappedTorus((pos-vec3( 0.0,0.30, 1.0))*vec3(1,-1,1), vec2(0.866025,-0.5), 0.25, 0.05), 25.0) );
    res = opU( res, vec2( sdBoxFrame(    pos-vec3( 0.0,0.25, 0.0), vec3(0.3,0.25,0.2), 0.025 ), 16.9 ) );
        res = opU( res, vec2( sdCone(        pos-vec3( 0.0,0.45,-1.0), vec2(0.6,0.8),0.45 ), 55.0 ) );
    res = opU( res, vec2( sdCappedCone(  pos-vec3( 0.0,0.25,-2.0), 0.25, 0.25, 0.1 ), 13.67 ) );
    res = opU( res, vec2( sdSolidAngle(  pos-vec3( 0.0,0.00,-3.0), vec2(3,4)/5.0, 0.4 ), 49.13 ) );
    }

    // bounding box
    if( sdBox( pos-vec3(1.0,0.3,-1.0),vec3(0.35,0.3,2.5) )<res.x )
    {
        res = opU( res, vec2( sdTorus(      (pos-vec3( 1.0,0.30, 1.0)).xzy, vec2(0.25,0.05) ), 7.1 ) );
    res = opU( res, vec2( sdBox(         pos-vec3( 1.0,0.25, 0.0), vec3(0.3,0.25,0.1) ), 3.0 ) );
    res = opU( res, vec2( sdCapsule(     pos-vec3( 1.0,0.00,-1.0),vec3(-0.1,0.1,-0.1), vec3(0.2,0.4,0.2), 0.1  ), 31.9 ) );
        res = opU( res, vec2( sdCylinder(    pos-vec3( 1.0,0.25,-2.0), vec2(0.15,0.25) ), 8.0 ) );
    res = opU( res, vec2( sdHexPrism(    pos-vec3( 1.0,0.2,-3.0), vec2(0.2,0.05) ), 18.4 ) );
    }

    // bounding box
    if( sdBox( pos-vec3(-1.0,0.35,-1.0),vec3(0.35,0.35,2.5))<res.x )
    {
        res = opU( res, vec2( sdPyramid(    pos-vec3(-1.0,-0.6,-3.0), 1.0 ), 13.56 ) );
        res = opU( res, vec2( sdOctahedron( pos-vec3(-1.0,0.15,-2.0), 0.35 ), 23.56 ) );
    res = opU( res, vec2( sdTriPrism(   pos-vec3(-1.0,0.15,-1.0), vec2(0.3,0.05) ),43.5 ) );
    res = opU( res, vec2( sdEllipsoid(  pos-vec3(-1.0,0.25, 0.0), vec3(0.2, 0.25, 0.05) ), 43.17 ) );
    res = opU( res, vec2( sdHorseshoe(  pos-vec3(-1.0,0.25, 1.0), vec2(cos(1.3),sin(1.3)), 0.2, 0.3, vec2(0.03,0.08) ), 11.5 ) );
    }

    // bounding box
    if( sdBox( pos-vec3(2.0,0.3,-1.0),vec3(0.35,0.3,2.5) )<res.x )
    {
    res = opU( res, vec2( sdOctogonPrism(pos-vec3( 2.0,0.2,-3.0), 0.2, 0.05), 51.8 ) );
    res = opU( res, vec2( sdCylinder(    pos-vec3( 2.0,0.14,-2.0), vec3(0.1,-0.1,0.0), vec3(-0.2,0.35,0.1), 0.08), 31.2 ) );
        res = opU( res, vec2( sdCappedCone(  pos-vec3( 2.0,0.09,-1.0), vec3(0.1,0.0,0.0), vec3(-0.2,0.40,0.1), 0.15, 0.05), 46.1 ) );
    res = opU( res, vec2( sdRoundCone(   pos-vec3( 2.0,0.15, 0.0), vec3(0.1,0.0,0.0), vec3(-0.1,0.35,0.1), 0.15, 0.05), 51.7 ) );
    res = opU( res, vec2( sdRoundCone(   pos-vec3( 2.0,0.20, 1.0), 0.2, 0.1, 0.3 ), 37.0 ) );
    }

    return res;*/
}

// https://iquilezles.org/articles/boxfunctions
vec2 iBox( in vec3 ro, in vec3 rd, in vec3 rad )
{
    vec3 m = 1.0/rd;
    vec3 n = m*ro;
    vec3 k = abs(m)*rad;
    vec3 t1 = -n - k;
    vec3 t2 = -n + k;
        return vec2( max( max( t1.x, t1.y ), t1.z ),
                     min( min( t2.x, t2.y ), t2.z ) );
}

vec3 raycast( in vec3 ro, in vec3 rd , out vec3 normal )
{
    vec3 res = vec3(-1.0, -1.0, -1.0);

    float tmin = 1.0;
    float tmax = 20000.0;

    // Voxels
    vec3 texSize = textureSize(dataFieldTex, 0);

//    dist = sdVoxels(pos, rd, dataFieldTex);
//    vec4 hittingVoxel = getHittingVoxel(ro, rd, dataFieldTex);
    int invisibleMatIndex = 0; // 0 = air
    vec4 hittingVoxel = getHittingStackedMaterial(ro, rd, matIndicesTex, matHeightsTex, invisibleMatIndex, normal);
    if (hittingVoxel.w > -1) {
        tmax = length(hittingVoxel.xyz - ro);
        res = vec3(tmax, 5.0, -1);
//        res = vec3(tmax, 5.0, hittingVoxel.w);


//        vec3 normals[6];
//        normals[0] = vec3( 1., 0., 0.); normals[1] = vec3( -1.,  0.,  0.);
//        normals[2] = vec3( 0., 1., 0.); normals[3] = vec3(  0., -1.,  0.);
//        normals[4] = vec3( 0., 0., 1.); normals[5] = vec3(  0.,  0., -1.);
//        nor = normals[int(hittingVoxel.w)];
    }

    // raytrace floor plane
    float tp1 = (0.0-ro.y)/rd.y;
//    if( tp1>0.0 )
//    {
//        tmax = min( tmax, tp1 );
//        res = vec2( tp1, 1.0 );
//    }
    //else return res;

    // raymarch primitives
    vec2 tb = iBox( ro-vec3(0.0,0.4,-10.5)*10000.0, rd, vec3(2.5,0.41,3.0)*10000.0 );
//    if( tb.x<tb.y && tb.y>0.0 && tb.x<tmax)
//    {
        //return vec2(tb.x,2.0);
//        tmin = max(tb.x,tmin);
//        tmax = min(tb.y,tmax);

        float t = tmin;
        for( int i=0; i<70 && t<tmax; i++ )
        {
            hit_sphere hit;
            vec2 h = map( ro+rd*t , rd, hit);
            if( abs(h.x) < (0.0001*t) )
            {
                res = vec3(t, h.y, -1);
                normal = (h.y<1.5) ? vec3(0.0,1.0,0.0) : calcNormal( ro+rd*t );
                break;
            }
            t += h.x;
        }
//    }

    return res;
}

// https://iquilezles.org/articles/rmshadows
float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax )
{
    // bounding volume
    float tp = (0.8-ro.y)/rd.y; if( tp>0.0 ) tmax = min( tmax, tp );

    float res = 1.0;
    float t = mint;
    for( int i=0; i<24; i++ )
    {
        hit_sphere hit;
        float h = map( ro + rd*t , rd, hit).x;
        float s = clamp(8.0*h/t,0.0,1.0);
        res = min( res, s );
        t += clamp( h, 0.01, 0.2 );
        if( res<0.004 || t>tmax ) break;
    }
    res = clamp( res, 0.0, 1.0 );
    return res*res*(3.0-2.0*res);
}

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal( in vec3 pos )
{
    hit_sphere hit;
#if 0
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy*map( pos + e.xyy , e.xyy, hit).x +
                                          e.yyx*map( pos + e.yyx , e.yyx, hit).x +
                                          e.yxy*map( pos + e.yxy , e.yxy, hit).x +
                                          e.xxx*map( pos + e.xxx , e.xxx, hit).x );
#else
    // inspired by tdhooper and klems - a way to prevent the compiler from inlining map() 4 times
    vec3 n = vec3(0.0);
    for( int i=0; i<4; i++ )
    {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*map(pos+0.0005*e, e, hit).x;
      //if( n.x+n.y+n.z>100.0 ) break;
    }
    return normalize(n);
#endif
}

// https://iquilezles.org/articles/nvscene2008/rwwtt.pdf
float calcAO( in vec3 pos, in vec3 nor )
{
        float occ = 0.0;
    float sca = 1.0;
    hit_sphere hit;
    for( int i=0; i<5; i++ )
    {
        float h = 0.01 + 0.12*float(i)/4.0;
        float d = map( pos + h*nor , nor, hit).x;
        occ += (h-d)*sca;
        sca *= 0.95;
        if( occ>0.35 ) break;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 ) * (0.5+0.5*nor.y);
}

// https://iquilezles.org/articles/checkerfiltering
float checkersGradBox( in vec2 p, in vec2 dpdx, in vec2 dpdy )
{
    // filter kernel
    vec2 w = abs(dpdx)+abs(dpdy) + 0.001;
    // analytical integral (box filter)
    vec2 i = 2.0*(abs(fract((p-0.5*w)*0.5)-0.5)-abs(fract((p+0.5*w)*0.5)-0.5))/w;
    // xor pattern
    return 0.5 - 0.5*i.x*i.y;
}

vec3 render( in vec3 ro, in vec3 rd, in vec3 rdx, in vec3 rdy )
{
    // background
    vec3 col = vec3(0.7, 0.7, 0.9) - max(rd.y,0.0)*0.3;

    // raycast scene
    vec3 nor;
    vec3 res = raycast(ro,rd, nor);
    float t = res.x;
    float m = res.y;

    vec3 pos = ro + t*rd;
//    float stackEval = evaluateLayerStackImplicitSurface(pos, 5.0, 0, 1000);
//    return vec3(stackEval + 1.0) * .5;
//    return rd;

    if( m>-0.5 )
    {
        vec3 ref = reflect( rd, nor );

//        return abs(normalize(nor));

        // material
        col = 0.2 + 0.2*sin( m*2.0 + vec3(0.0,1.0,2.0) );
        float ks = 1.0;

        if( m<1.5 )
        {
            // project pixel footprint into the plane
            vec3 dpdx = ro.y*(rd/rd.y-rdx/rdx.y);
            vec3 dpdy = ro.y*(rd/rd.y-rdy/rdy.y);

            float f = checkersGradBox( 3.0*pos.xz, 3.0*dpdx.xz, 3.0*dpdy.xz );
            col = 0.15 + f*vec3(0.05);
            ks = 0.4;
        }

        // lighting
#if AMBIENT_OCCLUSION
        float occ = calcAO( pos, nor );
#else
        float occ = 1.0;
#endif

                vec3 lin = vec3(0.0);

        // sun
        {
            vec3  lig = normalize( vec3(-0.5, 0.4, -0.6) );
            vec3  hal = normalize( lig-rd );
            float dif = clamp( dot( nor, lig ), 0.0, 1.0 );
          //if( dif>0.0001 )
                      dif *= calcSoftshadow( pos, lig, 0.02, 2.5 );
                        float spe = pow( clamp( dot( nor, hal ), 0.0, 1.0 ),16.0);
                  spe *= dif;
                  spe *= 0.04+0.96*pow(clamp(1.0-dot(hal,lig),0.0,1.0),5.0);
                //spe *= 0.04+0.96*pow(clamp(1.0-sqrt(0.5*(1.0-dot(rd,lig))),0.0,1.0),5.0);
            lin += col*2.20*dif*vec3(1.30,1.00,0.70);
            lin +=     5.00*spe*vec3(1.30,1.00,0.70)*ks;
        }
        // sky
        {
            float dif = sqrt(clamp( 0.5+0.5*nor.y, 0.0, 1.0 ));
                  dif *= occ;
            float spe = smoothstep( -0.2, 0.2, ref.y );
                  spe *= dif;
                  spe *= 0.04+0.96*pow(clamp(1.0+dot(nor,rd),0.0,1.0), 5.0 );
          //if( spe>0.001 )
                  spe *= calcSoftshadow( pos, ref, 0.02, 2.5 );
            lin += col*0.60*dif*vec3(0.40,0.60,1.15);
            lin +=     2.00*spe*vec3(0.40,0.60,1.30)*ks;
        }
        // back
        {
                float dif = clamp( dot( nor, normalize(vec3(0.5,0.0,0.6))), 0.0, 1.0 )*clamp( 1.0-pos.y,0.0,1.0);
                  dif *= occ;
                lin += col*0.55*dif*vec3(0.25,0.25,0.25);
        }
        // sss
        {
            float dif = pow(clamp(1.0+dot(nor,rd),0.0,1.0),2.0);
                  dif *= occ;
                lin += col*0.25*dif*vec3(1.00,1.00,1.00);
        }

                col = lin;

//        col = mix( col, vec3(0.7,0.7,0.9), 1.0-exp( -0.0001*t*t*t ) );
    }

        return vec3( clamp(col,0.0,1.0) );
}

void main()
{
    float time = 0;
    vec3 mo = vec3(0.0);
    // camera
//    vec3 ta = vec3( 0.25, -0.75, -0.75 );
//    vec3 ro = ta + vec3( 4.5*cos(0.1*time + 7.0*mo.x), 2.2, 4.5*sin(0.1*time + 7.0*mo.x) );

    mat4 matrix = scene.proj_matrix * scene.mv_matrix;
    vec3 ro = vec4(matrix * vec4(0, 0, 0, 1)).xyz;
    vec3 ta = ro - rayDirection;
//    vec3 ta = ro + vec3(1, -1, -1);
    // camera-to-world transformation
//    mat3 ca = setCamera( ro, ta, 0.0 );
//    vec3 ro;
//    mat3 ca = setCamera(ro);
    mat4 ca = scene.proj_matrix * scene.mv_matrix;
//    ro = vec4(ca * vec4(0, 0, 0, 1)).xyz;
//    vec3 fwd =  vec4(ca * vec4(0, 0, 1, 1)).xyz - ro;
//    FragColor = vec4(ro/300.0, .80);
//    return;
    ro = scene.camera_pos.xyz;

    vec3 iResolution = vec3(scene.canvas_width, scene.canvas_height, 1.0);
    vec2 fragCoord = gl_FragCoord.xy;

    vec3 tot = vec3(0.0);

    for( int m=0; m<NB_SAMPLES; m++ ) {
        for( int n=0; n<NB_SAMPLES; n++ )
        {
            // pixel coordinates
            vec2 o = vec2(float(m),float(n)) / float(NB_SAMPLES) - 0.5;
            vec2 p = (2.0*(fragCoord+o)-iResolution.xy)/iResolution.y;

            // focal length
            const float fl = 2.5;

            // ray direction
//            vec3 rd = ca * normalize( vec3(p,fl) );
            vec3 rd = rayDirection;
//            vec3 rd = vec4(ca * vec4(normalize( vec3(p,fl) ), 1.0)).xyz;

             // ray differentials
            vec2 px = (2.0*(fragCoord+vec2(1.0,0.0))-iResolution.xy)/iResolution.y;
            vec2 py = (2.0*(fragCoord+vec2(0.0,1.0))-iResolution.xy)/iResolution.y;
//            vec3 rdx = ca * normalize( vec3(px,fl) );
//            vec3 rdy = ca * normalize( vec3(py,fl) );
            vec3 rdx = rd + vec3(0.01, 0, 0);
            vec3 rdy = rd + vec3(0, 0, 0.01);
//            vec3 rdx = vec4(ca * vec4(normalize( vec3(px,fl) ), 1.0)).xyz;
//            vec3 rdy = vec4(ca * vec4(normalize( vec3(py,fl) ), 1.0)).xyz;

            // render
            vec3 col = render( ro, rd, rdx, rdy );

            // gain
            // col = col*3.0/(2.5+col);

                    // gamma
            col = pow( col, vec3(0.4545) );

            tot += col;
        }
    }
    tot /= float(NB_SAMPLES*NB_SAMPLES);

    FragColor = vec4( tot, 1.0 );
}




