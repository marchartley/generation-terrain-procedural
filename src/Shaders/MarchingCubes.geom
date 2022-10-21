/**** Geometry Shader Marching Cubes
* Copyright Cyril Crassin, Junuary 2007.
* This code is partially based on the example of
* Paul Bourke \”Polygonising a scalar field\” located at :
* http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
****/

//GLSL version 1.20
#version 430
//New G80 extensions
#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable

layout ( points ) in;
layout ( triangle_strip, max_vertices = 16 ) out;

//Volume data field texture
uniform sampler3D dataFieldTex;
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

uniform bool useMarchingCubes;

//Vertices position for fragment shader
in vec3 initialVertPos[];
in vec3 realNormal[];

out vec3 ginitialVertPos;
out vec3 grealNormal;
out vec4 gcolor;

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
    if (pos.z <= 1)
        val = 1.0;
    val = (val < min_isolevel || val > max_isolevel) ? -1000 : val;
    return val - 0.5;
}
//Get vertex i value within current marching cube
float cubeVal(vec3 voxelPos, int i){
    return cubeVal(cubePos(voxelPos, i));
}

//Get triangle table value
int triTableValue(int i, int j){
    return texelFetch2D(triTableTex, ivec2(j, i), 0).a;
}

//Compute interpolated vertex along an edge
vec3 vertexInterp(float isolevel, vec3 v0, float l0, vec3 v1, float l1){
    return mix(v0, v1, (isolevel-l0)/(l1-l0));
}

int getCubeIndex(vec3 voxPos, out vec3 normal) {
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

    normal = vec3(0, 0, 0);

    if (cubeindex != 0 && cubeindex != 255) {
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


        vec3 edge1 = vertlist[triTableValue(cubeindex, 0)] - vertlist[triTableValue(cubeindex, 1)];
        vec3 edge2 = vertlist[triTableValue(cubeindex, 0)] - vertlist[triTableValue(cubeindex, 2)];
        normal = normalize(cross(edge1, edge2));
    }
    return cubeindex;
}

//Geometry Shader entry point
void main(void) {
    vec4 position = vec4(initialVertPos[0], 1.0);
    vec3 voxPos = position.xyz;
    gcolor = vec4(1, 1, 1, 1);

    gl_Position = proj_matrix * mv_matrix * position;

    if (!useMarchingCubes) {
        if (cubeVal(position.xyz) < isolevel) return;
        ginitialVertPos = position.xyz;
        bool nTop    = (cubeVal(position.xyz + vec3( 0,  0,  1)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 0,  0,  1)) < -1000);
        bool nBottom = (cubeVal(position.xyz + vec3( 0,  0, -1)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 0,  0, -1)) < -1000);
        bool nRight  = (cubeVal(position.xyz + vec3( 1,  0,  0)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 1,  0,  0)) < -1000);
        bool nLeft   = (cubeVal(position.xyz + vec3(-1,  0,  0)) < isolevel); // ^^ cubeVal(position.xyz + vec3(-1,  0,  0)) < -1000);
        bool nFront  = (cubeVal(position.xyz + vec3( 0, -1,  0)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 0, -1,  0)) < -1000);
        bool nBack   = (cubeVal(position.xyz + vec3( 0,  1,  0)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 0,  1,  0)) < -1000);
        // Front
        grealNormal = vec3(0, -1, 0);
        if (nFront) {
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 0, 0, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 0, 1, 0));
            EmitVertex();
        }
        // Left-front corner
        if (nFront || nRight) {
            gl_Position = proj_matrix * mv_matrix * (position + vec4(1, 0, 0, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(1, 0, 1, 0));
            EmitVertex();
        }
        // Right
        grealNormal = vec3(1, 0, 0);
        if (nRight || nBack) {
            gl_Position = proj_matrix * mv_matrix * (position + vec4(1, 1, 0, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(1, 1, 1, 0));
            EmitVertex();
        }
        // Right
        grealNormal = vec3(0, 1, 0);
        if (nBack || nLeft) {
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 1, 0, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 1, 1, 0));
            EmitVertex();
        }
        // Back
        grealNormal = vec3(-1, 0, 0);
        if (nLeft) {
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 0, 0, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 0, 1, 0));
            EmitVertex();
        }
        EndPrimitive();

        // Bottom
        grealNormal = vec3(0, 0, -1);
        if (nBottom) {
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 0, 0, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(1, 0, 0, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 1, 0, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(1, 1, 0, 0));
            EmitVertex();
            EndPrimitive();
        }

        // Top

        if (nTop) {
            grealNormal = vec3(0, 0, 1);
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 0, 1, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(1, 0, 1, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(0, 1, 1, 0));
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * (position + vec4(1, 1, 1, 0));
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
        int cubeindex = getCubeIndex(voxPos, normal);

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
                gl_Position = proj_matrix * mv_matrix * position;
                ginitialVertPos = position.xyz;
                grealNormal = normal;
                EmitVertex();

                position= vec4(vertlist[triTableValue(cubeindex, i+1)], 1);
                gl_Position = proj_matrix * mv_matrix * position;
                ginitialVertPos = position.xyz;
                grealNormal = normal;
                EmitVertex();

                position= vec4(vertlist[triTableValue(cubeindex, i+2)], 1);
                gl_Position = proj_matrix * mv_matrix * position;
                ginitialVertPos = position.xyz;
                grealNormal = normal;
                EmitVertex();
                EndPrimitive();
            }else{
                break;
            }

            i = i + 3;
        }
    }
}
