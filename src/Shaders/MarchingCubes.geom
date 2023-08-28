/**** Geometry Shader Marching Cubes
* Copyright Cyril Crassin, Junuary 2007.
* This code is partially based on the example of
* Paul Bourke \”Polygonising a scalar field\” located at :
* http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
****/

//GLSL version 1.20
#version 430
//New G80 extensions
//#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable

layout ( points ) in;
layout ( triangle_strip, max_vertices = 24 ) out;

//Volume data field texture
uniform sampler3D dataFieldTex;
uniform sampler3D dataChangesFieldTex;
//Edge table texture
uniform isampler2D edgeTableTex;
//Triangles table texture
uniform isampler2D triTableTex;

//uniform sampler2D allBiomesColorTextures;

//Global iso level
uniform float isolevel;
//Marching cubes vertices decal
uniform vec3 vertDecals[8];

uniform float min_isolevel = -1000.0;
uniform float max_isolevel =  1000.0;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec3 min_vertice_positions;
uniform vec3 max_vertice_positions;

uniform float offsetX;
uniform float offsetY;
uniform float offsetZ;
uniform vec3 scale = vec3(1.0);
uniform vec3 rotation = vec3(0.0, 0.0, 0.0);

uniform bool useMarchingCubes;

//Vertices position for fragment shader
in vec3 initialVertPos[];
in vec3 realNormal[];

out vec3 ginitialVertPos;
out vec3 grealNormal;
out vec4 gcolor;
out float gdensity;
//out vec4 gPosition;

uniform int voxels_displayed_on_borders = 1;

vec4 rotate(vec4 pos, vec3 rot, vec3 center) {
    return pos;
    mat3 Rx = mat3(
        1, 0, 0,
        0, cos(rot.x), sin(rot.x),
        0, -sin(rot.x), cos(rot.x)
        );
    mat3 Ry = mat3(
        cos(rot.y), 0, -sin(rot.y),
        0, 1, 0,
        sin(rot.y), 0, cos(rot.y)
        );
    mat3 Rz = mat3(
        cos(rot.z), sin(rot.z), 0,
        -sin(rot.z), cos(rot.z), 0,
        0, 0, 1
        );
    return vec4(Rx * Ry * Rz * (pos.xyz - center) + center, 1);
}

//Get vertex i position within current marching cube
vec3 cubePos(vec3 voxelPos, int i){
    return vec4(vec4(voxelPos + vertDecals[i], 1.0)).xyz;
}

//Get vertex i value within current marching cube
float cubeVal(vec3 pos){
    vec3 texSize = textureSize(dataFieldTex, 0);
    if (pos.x == 0 || pos.y == 0 || pos.z == 0) return -1;
    if (pos.x < min_vertice_positions.x || max_vertice_positions.x < pos.x ||
        pos.y < min_vertice_positions.y || max_vertice_positions.y < pos.y ||
        pos.z < min_vertice_positions.z || max_vertice_positions.z < pos.z) return -1;
    if (pos.x >= texSize.x || pos.y >= texSize.y || pos.z >= texSize.z) return -1;
    float val = texture(dataFieldTex, pos/texSize).a;
//    if (pos.z <= 1)
//        val = 1.0;
    val = (val < min_isolevel || val > max_isolevel) ? -1 : val;
    return val - 0.5;
}
//Get vertex i value within current marching cube
float cubeVal(vec3 voxelPos, int i){
    return cubeVal(cubePos(voxelPos, i));
}

vec3 getNormal(vec3 p) {
    vec3 maxDims = min(textureSize(dataFieldTex, 0), max_vertice_positions);
    vec3 minDims = max(vec3(0.0), min_vertice_positions);
    vec3 x0 = clamp(p + vec3(0.1, 0, 0), minDims, maxDims);
    vec3 x1 = clamp(p - vec3(0.1, 0, 0), minDims, maxDims);
    vec3 y0 = clamp(p + vec3(0, 0.1, 0), minDims, maxDims);
    vec3 y1 = clamp(p - vec3(0, 0.1, 0), minDims, maxDims);
    vec3 z0 = clamp(p + vec3(0, 0, 0.1), minDims, maxDims);
    vec3 z1 = clamp(p - vec3(0, 0, 0.1), minDims, maxDims);
    float dx = cubeVal(x0) - cubeVal(x1) / length(x1 - x0);
    float dy = cubeVal(y0) - cubeVal(y1) / length(y1 - y0);
    float dz = cubeVal(z0) - cubeVal(z1) / length(z1 - z0);
//    dx = clamp(dx, -1, 1);
//    dy = clamp(dy, -1, 1);
//    dz = clamp(dz, -1, 1);
    /*
    float dx = pow(clamp(cubeVal(x0) - cubeVal(x1), -1, 1), 5);
    float dy = pow(clamp(cubeVal(y0) - cubeVal(y1), -1, 1), 5);
    float dz = pow(clamp(cubeVal(z0) - cubeVal(z1), -1, 1), 5);
    */
    return normalize(-vec3(dx, dy, dz));
}

//Get triangle table value
int triTableValue(int i, int j){
    return texelFetch2D(triTableTex, ivec2(j, i), 0).a;
}

//Compute interpolated vertex along an edge
vec3 vertexInterp(float isolevel, vec3 v0, float l0, vec3 v1, float l1){
    return mix(v0, v1, clamp((isolevel-l0)/(l1-l0), 0.0, 1.0)) * scale + vec3(offsetX, offsetY, offsetZ);
}

float mincomp(vec2 v) { return min(v.x, v.y); }
float maxcomp(vec2 v) { return max(v.x, v.y); }
float mincomp(vec3 v) { return min(min(v.x, v.y), v.z); }
float maxcomp(vec3 v) { return max(max(v.x, v.y), v.z); }
vec4 getPosition(vec4 position, vec3 _offset, vec3 center)
{
//    return position + vec4(_offset, 0.0);
    _offset += vec3(offsetX, offsetY, offsetZ);
    position *= vec4(scale, 1.0);

    float distToLimits = (voxels_displayed_on_borders > 1 ? min(mincomp(abs(position.xyz - min_vertice_positions)), mincomp(abs(position.xyz + vec3(1.0) - max_vertice_positions))) : 1.0);
    vec3 off = _offset * (clamp(distToLimits / float(voxels_displayed_on_borders), 0.0, 1.0));
    return rotate(clamp(position + vec4(off, 0.0), vec4(min_vertice_positions, 1.0), vec4(max_vertice_positions, 1.0)), rotation, center);
//    return clamp (position + vec4(_offset, 0.0), vec4(min_vertice_positions, 1.0), vec4(max_vertice_positions, 1.0));
}

int getCubeIndex(vec3 voxPos) {
    int cubeindex = 0;
    float cubeVal0 = cubeVal(voxPos, 0);
    float cubeVal1 = cubeVal(voxPos, 1);
    float cubeVal2 = cubeVal(voxPos, 2);
    float cubeVal3 = cubeVal(voxPos, 3);
    float cubeVal4 = cubeVal(voxPos, 4);
    float cubeVal5 = cubeVal(voxPos, 5);
    float cubeVal6 = cubeVal(voxPos, 6);
    float cubeVal7 = cubeVal(voxPos, 7);
    float refined_isolevel = isolevel + 0.0001;
    //Determine the index into the edge table which
    //tells us which vertices are inside of the surface
    cubeindex  = int(cubeVal0 < refined_isolevel);
    cubeindex += int(cubeVal1 < refined_isolevel)*2;
    cubeindex += int(cubeVal2 < refined_isolevel)*4;
    cubeindex += int(cubeVal3 < refined_isolevel)*8;
    cubeindex += int(cubeVal4 < refined_isolevel)*16;
    cubeindex += int(cubeVal5 < refined_isolevel)*32;
    cubeindex += int(cubeVal6 < refined_isolevel)*64;
    cubeindex += int(cubeVal7 < refined_isolevel)*128;
    return cubeindex;
}


float getDensity(vec3 pos/*, float resolution*/) {
    vec3 texSize = vec3(textureSize(dataFieldTex, 0));
    vec3 offsets = vec3(0.0); // 0.5 / texSize;
    float density = 0; // texture(dataFieldTex, pos).a;
    float surrounding = 2.f;
    for (float x = 0; x < surrounding + 1.0; x += 1.0) {
        for (float y = 0; y < surrounding + 1.0; y += 1.0) {
            for (float z = 0; z < surrounding + 1.0; z += 1.0) {
                float divisor = 0.1; //float(int(surrounding*.5f));
                vec3 newPos = pos + vec3(x - surrounding*.5f, y - surrounding*.5f, z - surrounding*.5f) / texSize; // * (resolution / divisor);
//                vec3 newPos = pos + vec3(resolution * (x/surrounding - surrounding * .5f), resolution * (y/surrounding - surrounding * .5f), resolution * (z/surrounding - surrounding * .5f));
                float val = texture(dataFieldTex, newPos - offsets).a;
                density = max(density, val);
            }
        }
    }
    return (density - 0.5);
}

//Geometry Shader entry point
void main(void) {
    vec4 position = vec4(initialVertPos[0], 1.0);
    vec3 voxPos = position.xyz;
    vec3 dataSize = vec3(textureSize(dataFieldTex, 0));
    vec3 center = dataSize * 0.5;

    vec3 changeSize = vec3(textureSize(dataChangesFieldTex, 0));
    vec3 evalPos = (changeSize.z > 1 ? voxPos / changeSize : vec3(voxPos.xy, 0.5) / changeSize);
    float changeVal = texture(dataChangesFieldTex, evalPos).a - 2.0;
    if (changeVal < -.2)
        gdensity = 10.0;
    else if (changeVal > .2)
        gdensity = 0.1999;
    else
        gdensity = 0.65; // getDensity(voxPos);

    if (!useMarchingCubes) {
        if (cubeVal(position.xyz) < isolevel)
            return;
        if (position.x + 1 < min_vertice_positions.x || position.x > max_vertice_positions.x ||
            position.y + 1 < min_vertice_positions.y || position.y > max_vertice_positions.y ||
            position.z + 1 < min_vertice_positions.z || position.z > max_vertice_positions.z)
            return;
        bool nTop    = (cubeVal(position.xyz + vec3( 0,  0,  1)) < isolevel) || position.z + (voxels_displayed_on_borders + 1) > max_vertice_positions.z;
        bool nBottom = (cubeVal(position.xyz + vec3( 0,  0, -1)) < isolevel) || position.z - (voxels_displayed_on_borders + 1) < min_vertice_positions.z;
        bool nRight  = (cubeVal(position.xyz + vec3( 1,  0,  0)) < isolevel) || position.x + (voxels_displayed_on_borders + 1) > max_vertice_positions.x;
        bool nLeft   = (cubeVal(position.xyz + vec3(-1,  0,  0)) < isolevel) || position.x - (voxels_displayed_on_borders + 1) < min_vertice_positions.x;
        bool nFront  = (cubeVal(position.xyz + vec3( 0, -1,  0)) < isolevel) || position.y - (voxels_displayed_on_borders + 1) < min_vertice_positions.y;
        bool nBack   = (cubeVal(position.xyz + vec3( 0,  1,  0)) < isolevel) || position.y + (voxels_displayed_on_borders + 1) > max_vertice_positions.y;

        // Front
        grealNormal = vec3(0, -1, 0);
        if (nFront) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 1), center).xyz;
//            gdensity =
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 0), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 0), center).xyz;
            EmitVertex();
            EndPrimitive();
        }
        // Right
        grealNormal = vec3(1, 0, 0);
        if (nRight) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 0), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 0), center).xyz;
            EmitVertex();
            EndPrimitive();
        }
        // Back
        grealNormal = vec3(0, 1, 0);
        if (nBack) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 0), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 0), center).xyz;
            EmitVertex();
            EndPrimitive();
        }
        // Left
        grealNormal = vec3(-1, 0, 0);
        if (nLeft) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 0), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 0), center).xyz;
            EmitVertex();
            EndPrimitive();
        }

        // Bottom
        grealNormal = vec3(0, 0, -1);
        if (nBottom) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 0), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 0), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 0), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 0), center).xyz;
            EmitVertex();
            EndPrimitive();
        }

        // Top
        if (nTop) {
            grealNormal = vec3(0, 0, 1);
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 1), center).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 1), center).xyz;
            EmitVertex();
            EndPrimitive();
        }
    } else {
        float cubeVal0 = cubeVal(voxPos, 0);
        float cubeVal1 = cubeVal(voxPos, 1);
        float cubeVal2 = cubeVal(voxPos, 2);
        float cubeVal3 = cubeVal(voxPos, 3);
        float cubeVal4 = cubeVal(voxPos, 4);
        float cubeVal5 = cubeVal(voxPos, 5);
        float cubeVal6 = cubeVal(voxPos, 6);
        float cubeVal7 = cubeVal(voxPos, 7);
        float refined_isolevel = isolevel + 0.0001;

        vec3 normal;
        int cubeindex = getCubeIndex(voxPos);

        //Cube is entirely in/out of the surface
        if (cubeindex == 0 || cubeindex == 255)
            return;

        vec3 vertlist[12];
        //Find the vertices where the surface intersects the cube
        vertlist[0] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 1), cubeVal1);
        vertlist[1] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 2), cubeVal2);
        vertlist[2] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 3), cubeVal3);
        vertlist[3] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 0), cubeVal0);
        vertlist[4] = vertexInterp(refined_isolevel, cubePos(voxPos, 4), cubeVal4, cubePos(voxPos, 5), cubeVal5);
        vertlist[5] = vertexInterp(refined_isolevel, cubePos(voxPos, 5), cubeVal5, cubePos(voxPos, 6), cubeVal6);
        vertlist[6] = vertexInterp(refined_isolevel, cubePos(voxPos, 6), cubeVal6, cubePos(voxPos, 7), cubeVal7);
        vertlist[7] = vertexInterp(refined_isolevel, cubePos(voxPos, 7), cubeVal7, cubePos(voxPos, 4), cubeVal4);
        vertlist[8] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 4), cubeVal4);
        vertlist[9] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 5), cubeVal5);
        vertlist[10] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 6), cubeVal6);
        vertlist[11] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 7), cubeVal7);

        int i = 0;
        while(true){
            if(triTableValue(cubeindex, i) != -1){
                position= vec4(vertlist[triTableValue(cubeindex, i)], 1);
                gl_Position = proj_matrix * mv_matrix * position; // getPosition(position, vec3(0), center);
                ginitialVertPos = position.xyz;
                grealNormal = getNormal(position.xyz);
                EmitVertex();

                position= vec4(vertlist[triTableValue(cubeindex, i+1)], 1);
                gl_Position = proj_matrix * mv_matrix * position; // getPosition(position, vec3(0), center);
                ginitialVertPos = position.xyz;
                grealNormal = getNormal(position.xyz);
                EmitVertex();

                position= vec4(vertlist[triTableValue(cubeindex, i+2)], 1);
                gl_Position = proj_matrix * mv_matrix *  position; // getPosition(position, vec3(0), center);
                ginitialVertPos = position.xyz;
                grealNormal = getNormal(position.xyz);
                EmitVertex();
                EndPrimitive();
            }else{
                break;
            }

            i = i + 3;
        }
    }
}
