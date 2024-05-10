#version 430
//New G80 extensions
//#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable

layout ( points ) in;
layout ( triangle_strip, max_vertices = 37) out;

//height data field texture
uniform sampler2D heightmapFieldTex;
uniform float maxHeight;
//uniform sampler2D allBiomesColorTextures;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec3 min_vertice_positions;
uniform vec3 max_vertice_positions;
uniform float ambiantOcclusionFactor = 0.0;
uniform float heightFactor = 0.1;
//Vertices position for fragment shader
in vec3 initialVertPos[];
in vec3 realNormal[];

out vec3 ginitialVertPos;
out vec3 grealNormal;
out vec4 gcolor;
out float gambiantOcclusion;


float getHeight(vec2 pos) {
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    return (pos.x >= texSize.x || pos.y >= texSize.y) ? 0 : texture(heightmapFieldTex, pos/texSize).a * maxHeight;
}
float getDisplacementLength(vec2 pos) {
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    return (pos.x >= texSize.x || pos.y >= texSize.y) ? 0 : texture(heightmapFieldTex, pos/texSize).b;
}

vec3 getNormal(vec2 pos) {
    vec3 normal = vec3(0, 0, 1);

    vec3 n0 = vec3(pos + vec2( 0,  0), getHeight(pos + vec2( 0,  0)));

    vec3 n1 = vec3(pos + vec2(-1, -1), getHeight(pos + vec2(-1, -1)));
    vec3 n2 = vec3(pos + vec2( 0, -1), getHeight(pos + vec2( 0, -1)));
    vec3 n3 = vec3(pos + vec2( 1, -1), getHeight(pos + vec2( 1, -1)));
    vec3 n4 = vec3(pos + vec2(-1,  0), getHeight(pos + vec2(-1,  0)));
    vec3 n5 = vec3(pos + vec2( 1,  0), getHeight(pos + vec2( 1,  0)));
    vec3 n6 = vec3(pos + vec2(-1,  1), getHeight(pos + vec2(-1,  1)));
    vec3 n7 = vec3(pos + vec2( 0,  1), getHeight(pos + vec2( 0,  1)));
    vec3 n8 = vec3(pos + vec2( 1,  1), getHeight(pos + vec2( 1,  1)));

    vec3 e1 = n1 - n0;
    vec3 e2 = n2 - n0;
    vec3 e3 = n3 - n0;
    vec3 e4 = n4 - n0;
    vec3 e5 = n5 - n0;
    vec3 e6 = n6 - n0;
    vec3 e7 = n7 - n0;
    vec3 e8 = n8 - n0;

    normal = cross(e1, e2) + cross(e2, e3) + cross(e3, e5) + cross(e5, e8) + cross(e8, e7) + cross(e7, e6) + cross(e6, e4) + cross(e4, e1);
    return normalize(normal);
}
float getAmbiantOcclusion(vec3 pos) {
    if (ambiantOcclusionFactor == 0.0) return 1.0;
    int nX = 20;
    int nY = 20;
    float fnX = float(nX);
    float fnY = float(nY);
    float pi = 3.141592;

    float occlusion = 0.0;
    float total = 0.0;
    for (int _r = 1; _r < 10; _r++) {
        for (int i = 0; i < nX; i++) {
            for (int j = 0; j < nY; j++) {
                float theta = (i / (fnX)) * 2.0 * pi;
                float phi = (j / (fnY - 1)) * pi - (pi * 0.5);
                float r = _r;

                vec3 ray = vec3(cos(theta) * cos(phi) * r, sin(theta) * cos(phi) * r, sin(phi) * r);
                occlusion += (pos.z + ray.z < getHeight((pos + ray).xy) ? 1.0 : 0.0);
                total += 1;
            }
        }
    }
    return pow(clamp(1.0 - smoothstep(0.0, 1.0, occlusion / total) + 0.5, 0.0, 1.0), 2.0);
}

void sendInfoVertex(vec4 vecPos) {
    float displacementStrength = 1.0;
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    grealNormal = getNormal(vecPos.xy);
    gcolor = vec4(texture(heightmapFieldTex, vecPos.xy/texSize).rgb, 1.0);
    vecPos.z = getHeight(vecPos.xy);
    vecPos += vec4(grealNormal * getDisplacementLength(vecPos.xy) * displacementStrength, 0.0);
    vecPos.z *= heightFactor;
    ginitialVertPos = vecPos.xyz;
    gl_Position = proj_matrix * mv_matrix * vecPos;
    EmitVertex();
}
void subdivision(vec2 vecPos) {
    vec4 v1 = vec4(vecPos + vec2(0, 0), getHeight(vecPos + vec2(0, 0)), 1.0);
    vec4 v2 = vec4(vecPos + vec2(1, 0), getHeight(vecPos + vec2(1, 0)), 1.0);
    vec4 v3 = vec4(vecPos + vec2(0, 1), getHeight(vecPos + vec2(0, 1)), 1.0);
    vec4 v4 = vec4(vecPos + vec2(1, 1), getHeight(vecPos + vec2(1, 1)), 1.0);
    vec4 v5 = vec4(vecPos + vec2(.5, .5), getHeight(vecPos + vec2(.5, .5)), 1.0);
/*
    sendInfoVertex(v1); sendInfoVertex(v2); sendInfoVertex(v5);
                        sendInfoVertex(v3); sendInfoVertex(v1);
    EndPrimitive();
    sendInfoVertex(v5); sendInfoVertex(v2); sendInfoVertex(v4);
                        sendInfoVertex(v3); sendInfoVertex(v5);
    EndPrimitive();*/

    sendInfoVertex(v1); sendInfoVertex(v2); sendInfoVertex(v5); sendInfoVertex(v4);
    EndPrimitive();
    sendInfoVertex(v4); sendInfoVertex(v3); sendInfoVertex(v5); sendInfoVertex(v1);
    EndPrimitive();
}
//Geometry Shader entry point
void main(void) {
    gambiantOcclusion = 1.0;
    float displacementStrength = 1.0;
    vec2 texSize = textureSize(heightmapFieldTex, 0);

    vec2 vecPos = initialVertPos[0].xy;
    if (vecPos.x >= texSize.x - 1 || vecPos.y >= texSize.y - 1) return;
    if (vecPos.x < min_vertice_positions.x || vecPos.x > max_vertice_positions.x ||
        vecPos.y < min_vertice_positions.y || vecPos.y > max_vertice_positions.y) return;

//    subdivision(vecPos);
//    return;


    /*vec4 v1 = vec4(vecPos + vec2(0, 0), getHeight(vecPos + vec2(0, 0)), 1.0);
    vec4 v2 = vec4(vecPos + vec2(1, 0), getHeight(vecPos + vec2(1, 0)), 1.0);
    vec4 v3 = vec4(vecPos + vec2(0, 1), getHeight(vecPos + vec2(0, 1)), 1.0);
    vec4 v4 = vec4(vecPos + vec2(1, 1), getHeight(vecPos + vec2(1, 1)), 1.0);

    sendInfoVertex(v1);
    sendInfoVertex(v2);
    sendInfoVertex(v3);
    sendInfoVertex(v4);
    */

    vec4 v11 = vec4(vecPos + vec2(0.00, 0.00), 0, 1.0);
    vec4 v12 = vec4(vecPos + vec2(0.25, 0.00), 0, 1.0);
    vec4 v13 = vec4(vecPos + vec2(0.50, 0.00), 0, 1.0);
    vec4 v14 = vec4(vecPos + vec2(0.75, 0.00), 0, 1.0);
    vec4 v15 = vec4(vecPos + vec2(1.00, 0.00), 0, 1.0);

    vec4 v21 = vec4(vecPos + vec2(0.00, 0.25), 0, 1.0);
    vec4 v22 = vec4(vecPos + vec2(0.25, 0.25), 0, 1.0);
    vec4 v23 = vec4(vecPos + vec2(0.50, 0.25), 0, 1.0);
    vec4 v24 = vec4(vecPos + vec2(0.75, 0.25), 0, 1.0);
    vec4 v25 = vec4(vecPos + vec2(1.00, 0.25), 0, 1.0);

    vec4 v31 = vec4(vecPos + vec2(0.00, 0.50), 0, 1.0);
    vec4 v32 = vec4(vecPos + vec2(0.25, 0.50), 0, 1.0);
    vec4 v33 = vec4(vecPos + vec2(0.50, 0.50), 0, 1.0);
    vec4 v34 = vec4(vecPos + vec2(0.75, 0.50), 0, 1.0);
    vec4 v35 = vec4(vecPos + vec2(1.00, 0.50), 0, 1.0);

    vec4 v41 = vec4(vecPos + vec2(0.00, 0.75), 0, 1.0);
    vec4 v42 = vec4(vecPos + vec2(0.25, 0.75), 0, 1.0);
    vec4 v43 = vec4(vecPos + vec2(0.50, 0.75), 0, 1.0);
    vec4 v44 = vec4(vecPos + vec2(0.75, 0.75), 0, 1.0);
    vec4 v45 = vec4(vecPos + vec2(1.00, 0.75), 0, 1.0);

    vec4 v51 = vec4(vecPos + vec2(0.00, 1.00), 0, 1.0);
    vec4 v52 = vec4(vecPos + vec2(0.25, 1.00), 0, 1.0);
    vec4 v53 = vec4(vecPos + vec2(0.50, 1.00), 0, 1.0);
    vec4 v54 = vec4(vecPos + vec2(0.75, 1.00), 0, 1.0);
    vec4 v55 = vec4(vecPos + vec2(1.00, 1.00), 0, 1.0);

    sendInfoVertex(v11);
    sendInfoVertex(v12);
    sendInfoVertex(v21);
    sendInfoVertex(v22);
    sendInfoVertex(v31);
    sendInfoVertex(v32);
    sendInfoVertex(v41);
    sendInfoVertex(v42);
    sendInfoVertex(v51);
    sendInfoVertex(v52);

    sendInfoVertex(v53);
    sendInfoVertex(v42);
    sendInfoVertex(v43);
    sendInfoVertex(v32);
    sendInfoVertex(v33);
    sendInfoVertex(v22);
    sendInfoVertex(v23);
    sendInfoVertex(v12);
    sendInfoVertex(v13);

    sendInfoVertex(v14);
    sendInfoVertex(v23);
    sendInfoVertex(v24);
    sendInfoVertex(v33);
    sendInfoVertex(v34);
    sendInfoVertex(v43);
    sendInfoVertex(v44);
    sendInfoVertex(v53);
    sendInfoVertex(v54);

    sendInfoVertex(v55);
    sendInfoVertex(v44);
    sendInfoVertex(v45);
    sendInfoVertex(v34);
    sendInfoVertex(v35);
    sendInfoVertex(v24);
    sendInfoVertex(v25);
    sendInfoVertex(v14);
    sendInfoVertex(v15);

    /*grealNormal = getNormal(v1.xy);
    gambiantOcclusion = getAmbiantOcclusion(v1.xyz);
    v1 += vec4(grealNormal * getDisplacementLength(v1.xy) * displacementStrength, 0.0);
    ginitialVertPos = v1.xyz;
    gcolor = vec4(texture(heightmapFieldTex, v1.xy/texSize).rgb, 1.0);
    v1.z *= heightFactor;
    gl_Position = proj_matrix * mv_matrix * v1;
    EmitVertex();

    grealNormal = getNormal(v2.xy);
    gambiantOcclusion = getAmbiantOcclusion(v2.xyz);
    v2 += vec4(grealNormal * getDisplacementLength(v2.xy) * displacementStrength, 0.0);
    ginitialVertPos = v2.xyz;
    gcolor = vec4(texture(heightmapFieldTex, v2.xy/texSize).rgb, 1.0);
    v2.z *= heightFactor;
    gl_Position = proj_matrix * mv_matrix * v2;
    EmitVertex();

    grealNormal = getNormal(v3.xy);
    gambiantOcclusion = getAmbiantOcclusion(v3.xyz);
    v3 += vec4(grealNormal * getDisplacementLength(v3.xy) * displacementStrength, 0.0);
    ginitialVertPos = v3.xyz;
    gcolor = vec4(texture(heightmapFieldTex, v3.xy/texSize).rgb, 1.0);
    v3.z *= heightFactor;
    gl_Position = proj_matrix * mv_matrix * v3;
    EmitVertex();

    grealNormal = getNormal(v4.xy);
    gambiantOcclusion = getAmbiantOcclusion(v4.xyz);
    v4 += vec4(grealNormal * getDisplacementLength(v4.xy) * displacementStrength, 0.0);
    ginitialVertPos = v4.xyz;
    gcolor = vec4(texture(heightmapFieldTex, v4.xy/texSize).rgb, 1.0);
    v4.z *= heightFactor;
    gl_Position = proj_matrix * mv_matrix * v4;
    EmitVertex();
    */
    EndPrimitive();
}
