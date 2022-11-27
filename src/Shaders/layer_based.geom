#version 430
//New G80 extensions
#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable

layout ( points ) in;
layout ( triangle_strip, max_vertices = 24 ) out; // Outputs a box

uniform sampler3D matIndicesTex;
uniform sampler3D matHeightsTex;

uniform float min_isolevel = -1000.0;
uniform float max_isolevel =  1000.0;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec3 min_vertice_positions;
uniform vec3 max_vertice_positions;

uniform bool useMarchingCubes;

uniform int voxels_displayed_on_borders = 5;

in vec3 initialVertPos[];

out vec3 ginitialVertPos;
out vec3 grealNormal;
out float gdensity;

float getDensityPerMaterial(int material)
{
    float density = 0.01;
    if (material == 0) {
        density = -10.0;
    } else if (material == 1) {
        density = -1.0;
    } else if (material == 2) {
        density = 0.25;
    } else if (material == 3) {
        density = 0.75;
    } else if (material == 4) {
        density = 1.4;
    } else if (material == 5) {
        density = 2.5;
    }
    return density;
}

bool isMatter(vec3 texPos)
{
    return texelFetch(matIndicesTex, ivec3(texPos), 0).a > 1;
}

float getMinHeight(vec3 texPos)
{
    ivec3 currentPos = ivec3(texPos.xy, texPos.z - 1);
    float currentHeight = 0.0;
    while (currentPos.z >= 0) {
        currentHeight += texelFetch(matHeightsTex, currentPos, 0).a;
        currentPos.z -= 1;
    }
    return currentHeight;
}

float getHeight(vec3 texPos)
{
    return texelFetch(matHeightsTex, ivec3(texPos), 0).a;
}

float mincomp(vec2 v) { return min(v.x, v.y); }
float maxcomp(vec2 v) { return max(v.x, v.y); }
float mincomp(vec3 v) { return min(min(v.x, v.y), v.z); }
float maxcomp(vec3 v) { return max(max(v.x, v.y), v.z); }
vec4 getPosition(vec4 position, vec3 _offset)
{
    float distToLimits = (voxels_displayed_on_borders > 1 ? min(mincomp(abs(position.xy - min_vertice_positions.xy)), mincomp(abs(position.xy + vec2(1.0) - max_vertice_positions.xy))) : 2.0);
    float factor = clamp(distToLimits / float(voxels_displayed_on_borders), 0.0, 2.0);
    vec3 off = _offset * vec3(factor, factor, 2.0) * .5 - vec3(.5, .5, 0.0);
    return position + vec4(.5, .5, .0, .0) + vec4(off, 0.0);
//    return clamp(position + vec4(.5, .5, .0, .0) + vec4(off, 0.0), vec4(min_vertice_positions, 1.0), vec4(max_vertice_positions, 1.0));
}

//Geometry Shader entry point
void main(void) {
    vec4 texture_position = vec4(initialVertPos[0].xyz, 1.0);
    vec4 position = vec4(texture_position.xy, 0.0, 1.0);

    gdensity = getDensityPerMaterial(int(round(texelFetch(matIndicesTex, ivec3(texture_position.xyz), 0).a)));

    if (gdensity <= 0.0) return; // Don't send to the frag shader if it's water/air
//    gdensity = max(0.0, gdensity + 0.5);

//    if (!useMarchingCubes) {
    if (true) {
        if (!isMatter(texture_position.xyz)) return;

        float height = getHeight(texture_position.xyz);
        if (height < 0.1) return;

        bool nTop    = true; // (cubeVal(position.xyz + vec3( 0,  0,  1)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 0,  0,  1)) < -1000);
        bool nBottom = true; // (cubeVal(position.xyz + vec3( 0,  0, -1)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 0,  0, -1)) < -1000);
        bool nRight  = true; // (cubeVal(position.xyz + vec3( 1,  0,  0)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 1,  0,  0)) < -1000);
        bool nLeft   = true; // (cubeVal(position.xyz + vec3(-1,  0,  0)) < isolevel); // ^^ cubeVal(position.xyz + vec3(-1,  0,  0)) < -1000);
        bool nFront  = true; // (cubeVal(position.xyz + vec3( 0, -1,  0)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 0, -1,  0)) < -1000);
        bool nBack   = true; // (cubeVal(position.xyz + vec3( 0,  1,  0)) < isolevel); // ^^ cubeVal(position.xyz + vec3( 0,  1,  0)) < -1000);

        float startHeight = getMinHeight(texture_position.xyz);
        float endHeight = startHeight + height;

        if (position.x + 1 < min_vertice_positions.x || position.x > max_vertice_positions.x ||
            position.y + 1 < min_vertice_positions.y || position.y > max_vertice_positions.y ||
            endHeight < min_vertice_positions.z || startHeight > max_vertice_positions.z)
            return;

        float minHeight = clamp(startHeight, min_vertice_positions.z, max_vertice_positions.z);
        float maxHeight = clamp(endHeight, min_vertice_positions.z, max_vertice_positions.z);

        // Front
        grealNormal = vec3(0, -1, 0);
        if (nFront) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, maxHeight));
            ginitialVertPos = getPosition(position, vec3(0, 0, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, minHeight));
            ginitialVertPos = getPosition(position, vec3(0, 0, minHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, maxHeight));
            ginitialVertPos = getPosition(position, vec3(1, 0, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, minHeight));
            ginitialVertPos = getPosition(position, vec3(1, 0, minHeight)).xyz;
            EmitVertex();
            EndPrimitive();
        }
        // Right
        grealNormal = vec3(1, 0, 0);
        if (nRight) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, maxHeight));
            ginitialVertPos = getPosition(position, vec3(1, 0, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, minHeight));
            ginitialVertPos = getPosition(position, vec3(1, 0, minHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, maxHeight));
            ginitialVertPos = getPosition(position, vec3(1, 1, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, minHeight));
            ginitialVertPos = getPosition(position, vec3(1, 1, minHeight)).xyz;
            EmitVertex();
            EndPrimitive();
        }
        // Right
        grealNormal = vec3(0, 1, 0);
        if (nBack) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, maxHeight));
            ginitialVertPos = getPosition(position, vec3(1, 1, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, minHeight));
            ginitialVertPos = getPosition(position, vec3(1, 1, minHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, maxHeight));
            ginitialVertPos = getPosition(position, vec3(0, 1, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, minHeight));
            ginitialVertPos = getPosition(position, vec3(0, 1, minHeight)).xyz;
            EmitVertex();
            EndPrimitive();
        }
        // Back
        grealNormal = vec3(-1, 0, 0);
        if (nLeft) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, maxHeight));
            ginitialVertPos = getPosition(position, vec3(0, 1, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, minHeight));
            ginitialVertPos = getPosition(position, vec3(0, 1, minHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, maxHeight));
            ginitialVertPos = getPosition(position, vec3(0, 0, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, minHeight));
            ginitialVertPos = getPosition(position, vec3(0, 0, minHeight)).xyz;
            EmitVertex();
            EndPrimitive();
        }

        // Bottom
        grealNormal = vec3(0, 0, -1);
        if (nBottom) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, minHeight));
            ginitialVertPos = getPosition(position, vec3(0, 0, minHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, minHeight));
            ginitialVertPos = getPosition(position, vec3(0, 1, minHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, minHeight));
            ginitialVertPos = getPosition(position, vec3(1, 0, minHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, minHeight));
            ginitialVertPos = getPosition(position, vec3(1, 1, minHeight)).xyz;
            EmitVertex();
            EndPrimitive();
        }

        // Top

        grealNormal = vec3(0, 0, 1);
        if (nTop) {
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 0, maxHeight));
            ginitialVertPos = getPosition(position, vec3(0, 0, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 0, maxHeight));
            ginitialVertPos = getPosition(position, vec3(1, 0, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(0, 1, maxHeight));
            ginitialVertPos = getPosition(position, vec3(0, 1, maxHeight)).xyz;
            EmitVertex();
            gl_Position = proj_matrix * mv_matrix * getPosition(position, vec3(1, 1, maxHeight));
            ginitialVertPos = getPosition(position, vec3(1, 1, maxHeight)).xyz;
            EmitVertex();
            EndPrimitive();
        }
    } /*else {
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
    }*/
}
