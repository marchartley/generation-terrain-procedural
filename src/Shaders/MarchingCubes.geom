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
uniform float ambiantOcclusionFactor = 0.0;
uniform float heightFactor = 1.0;

//Vertices position for fragment shader
in vec3 initialVertPos[];
in vec3 realNormal[];

out vec3 ginitialVertPos;
out vec3 grealNormal;
out vec4 gcolor;
out vec3 vColor;
out float gdensity;
out float gambiantOcclusion;
//out vec4 gPosition;

uniform int voxels_displayed_on_borders = 1;
uniform float waterRelativeHeight = 0.0;

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
    if (pos.x < 0 || pos.y < 0 || pos.z < 0) return -1;
    if (pos.x < min_vertice_positions.x || max_vertice_positions.x <= pos.x ||
        pos.y < min_vertice_positions.y || max_vertice_positions.y <= pos.y ||
        pos.z < min_vertice_positions.z || max_vertice_positions.z <= pos.z) return -1;
    if (pos.x >= texSize.x || pos.y >= texSize.y || pos.z >= texSize.z) return -1;

    // TEMP FIX
//    if (pos.z <= texSize.z * waterRelativeHeight) return 1.0;

    float val = texture(dataFieldTex, (pos + vec3(.5))/texSize).a;
//    if (pos.z <= 1)
//        val = 1.0;
    val = (val < min_isolevel || val > max_isolevel) ? -1 : val;
    return val - 0.5;
}
//Get vertex i value within current marching cube
float cubeVal(vec3 voxelPos, int i){
    return cubeVal(cubePos(voxelPos, i));
}

float __cubeVal(vec3 pos) {
    return cubeVal(pos);
}
vec3 getNormal(vec3 p) {
    if (p.x <= 1.0) return vec3(-1, 0, 0);
    else if (p.y <= 1.0) return vec3(0, -1, 0);
    else if (p.z <= 1.0) return vec3(0, 0, -1);

    else if (p.x >= max_vertice_positions.x - 1.0) return vec3(1, 0, 0);
    else if (p.y >= max_vertice_positions.y - 1.0) return vec3(0, 1, 0);
    else if (p.z >= max_vertice_positions.z - 1.0) return vec3(0, 0, 1);

    vec3 maxDims = min(textureSize(dataFieldTex, 0), max_vertice_positions);
    vec3 minDims = max(vec3(0.0), min_vertice_positions);

    // One central-difference stencil matches the CPU field gradient and keeps
    // fine features. The old two-radius average doubled all six texture reads.
    vec3 x0 = clamp(p + vec3(1,0,0), minDims, maxDims), x1 = clamp(p - vec3(1,0,0), minDims, maxDims);
    vec3 y0 = clamp(p + vec3(0,1,0), minDims, maxDims), y1 = clamp(p - vec3(0,1,0), minDims, maxDims);
    vec3 z0 = clamp(p + vec3(0,0,1), minDims, maxDims), z1 = clamp(p - vec3(0,0,1), minDims, maxDims);
    vec3 gradient = -vec3((cubeVal(x0) - cubeVal(x1)) / max(x0.x - x1.x, .000001),
                          (cubeVal(y0) - cubeVal(y1)) / max(y0.y - y1.y, .000001),
                          (cubeVal(z0) - cubeVal(z1)) / max(z0.z - z1.z, .000001));
    return dot(gradient, gradient) > 1e-8 ? normalize(gradient) : vec3(0,0,1);
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

int getCubeIndex(float samples[8]) {
    int index = 0;
    for (int corner = 0; corner < 8; ++corner)
        if (samples[corner] < isolevel + .0001) index |= 1 << corner;
    return index;
}

void setCellDensity(vec3 pos) {
    vec3 changeSize = vec3(textureSize(dataChangesFieldTex, 0));
    float changeVal = 0.0;
    if (all(greaterThan(changeSize, vec3(0)))) {
        vec3 evalPos = changeSize.z > 1.0 ? pos / changeSize : vec3(pos.xy, .5) / changeSize;
        changeVal = texture(dataChangesFieldTex, evalPos).a - 2.0;
    }
    gdensity = changeVal < -.2 ? 10.0 : (changeVal > .2 ? .1999 : .65);
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
                vec3 newPos = (pos + vec3(x - surrounding*.5f, y - surrounding*.5f, z - surrounding*.5f)) / texSize; // * (resolution / divisor);
//                vec3 newPos = pos + vec3(resolution * (x/surrounding - surrounding * .5f), resolution * (y/surrounding - surrounding * .5f), resolution * (z/surrounding - surrounding * .5f));
                float val = texture(dataFieldTex, newPos - offsets).a;
                density = max(density, val);
            }
        }
    }
    return (density - 0.5);
}

bool checkPos(vec3 pos) {
//    pos -= min_vertice_positions;
//    boxSize -= min_vertice_positions;
//    boxSize = clamp(boxSize, vec3(0, 0, 0), max_vertice_positions);
    return min_vertice_positions.x + 1 <= pos.x && pos.x + 1 <= max_vertice_positions.x &&
            min_vertice_positions.y + 1 <= pos.y && pos.y + 1 <= max_vertice_positions.y &&
            min_vertice_positions.z + 1 <= pos.z && pos.z + 1 <= max_vertice_positions.z;
}

float getAmbiantOcclusion(vec3 pos, vec3 normal) {
    if (ambiantOcclusionFactor == 0.0) return 1.0;
    vec3 texSize = vec3(textureSize(dataFieldTex, 0));
    vec3 offsets = vec3(0.0);
    float occlusion = 0.0;
    float total = 0.f;
    float surrounding = 4.f;
    vec3 up = normal;
    vec3 left = vec3(0.0, 1.0, 0.0);
    vec3 front = vec3(1.0, 0.0, 0.0);
    left = cross(up, front);
    if (!checkPos(pos - vec3(surrounding * .5, surrounding * .5, 0)) || !checkPos(pos + vec3(surrounding * .5, surrounding * .5, 0))) return 1.0;
    for (float x = -surrounding * .5f; x < surrounding * .5f + 1.0; x += 1.0) {
        for (float y = -surrounding * .5f; y < surrounding * .5f + 1.0; y += 1.0) {
            for (float z = -surrounding * .5f; z < surrounding * .5f + 1.0; z += 1.0) {
                vec3 off = x * front + y * left + z * up; //vec3(x - surrounding*.5f, y - surrounding*.5f, z - surrounding*.5f) * 1.0 + .5;
                vec3 newPos = (pos + off) / texSize;
                if (dot(off, normal) < 0) continue;
                occlusion += (texture(dataFieldTex, newPos - offsets).a > 0.5 && checkPos(newPos * texSize) ? 1.0 : 0.0);
                total += 1.0;
            }
        }
    }
    return  1.0 - occlusion / total;
}

//Geometry Shader entry point
void main(void) {
    vec4 position = vec4(initialVertPos[0], 1.0);
    vec3 voxPos = position.xyz;
    vec3 dataSize = vec3(textureSize(dataFieldTex, 0));
    vec3 center = dataSize * 0.5;


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
        setCellDensity(voxPos);
        grealNormal = vec3(0, -1, 0);
        if (nFront) {
            gambiantOcclusion = getAmbiantOcclusion(position.xyz + grealNormal * .5, grealNormal);
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 1), center).xyz;
//            gdensity =
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            EndPrimitive();
        }
        // Right
        grealNormal = vec3(1, 0, 0);
        if (nRight) {
            gambiantOcclusion = getAmbiantOcclusion(position.xyz + grealNormal * .5, grealNormal);
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            EndPrimitive();
        }
        // Back
        grealNormal = vec3(0, 1, 0);
        if (nBack) {
            gambiantOcclusion = getAmbiantOcclusion(position.xyz + grealNormal * .5, grealNormal);
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            EndPrimitive();
        }
        // Left
        if (nLeft) {
            grealNormal = vec3(-1, 0, 0);
            gambiantOcclusion = getAmbiantOcclusion(position.xyz + grealNormal * .5, grealNormal);
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            EndPrimitive();
        }

        // Bottom
        if (nBottom) {
            grealNormal = vec3(0, 0, -1);
            gambiantOcclusion = getAmbiantOcclusion(position.xyz + grealNormal * .5, grealNormal);
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 0), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 0), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            EndPrimitive();
        }

        // Top
        grealNormal = vec3(0, 0, -1);
        if (nTop) {
            gambiantOcclusion = getAmbiantOcclusion(position.xyz + grealNormal * .5, grealNormal);
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 0, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 0, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(0, 1, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, 1), center);
            ginitialVertPos = getPosition(position, vec3(1, 1, 1), center).xyz;
            vColor = vec3(0.0);
            EmitVertex();
            EndPrimitive();
        }
    } else {
        float samples[8];
        for (int corner = 0; corner < 8; ++corner) samples[corner] = cubeVal(voxPos, corner);
        int cubeindex = getCubeIndex(samples);
        if (cubeindex == 0 || cubeindex == 255) return;
        setCellDensity(voxPos);
        const ivec2 edgeCorners[12] = ivec2[12](
            ivec2(0,1), ivec2(1,2), ivec2(2,3), ivec2(3,0),
            ivec2(4,5), ivec2(5,6), ivec2(6,7), ivec2(7,4),
            ivec2(0,4), ivec2(1,5), ivec2(2,6), ivec2(3,7));
        vec3 positions[12];
        for (int edge = 0; edge < 12; ++edge) {
            ivec2 ends = edgeCorners[edge];
            positions[edge] = vertexInterp(isolevel + .0001,
                cubePos(voxPos, ends.x), samples[ends.x], cubePos(voxPos, ends.y), samples[ends.y]);
        }
        // Keep the geometry stage's live register set small. Large per-edge
        // caches reduce occupancy on mobile GPUs more than they save in reads.
        for (int triangle = 0; triangle < 15; triangle += 3) {
            int first = triTableValue(cubeindex, triangle);
            if (first < 0) break;
            for (int corner = 0; corner < 3; ++corner) {
                int edge = corner == 0 ? first : triTableValue(cubeindex, triangle + corner);
                ginitialVertPos = positions[edge];
                grealNormal = getNormal(ginitialVertPos);
                gambiantOcclusion = getAmbiantOcclusion(ginitialVertPos + grealNormal * .5, grealNormal);
                vec4 scaledPosition = vec4(ginitialVertPos, 1.0);
                scaledPosition.z *= heightFactor;
                gl_Position = proj_matrix * mv_matrix * scaledPosition;
                vColor = vec3(0.0);
                EmitVertex();
            }
            EndPrimitive();
        }
    }
}
