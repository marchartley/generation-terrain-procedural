#version 430
//New G80 extensions
//#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable

layout ( points ) in;
layout ( triangle_strip, max_vertices = 16) out;

//height data field texture
uniform sampler2D heightmapFieldTex;
uniform float maxHeight;
//uniform sampler2D allBiomesColorTextures;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec3 min_vertice_positions;
uniform vec3 max_vertice_positions;
//Vertices position for fragment shader
in vec3 initialVertPos[];
in vec3 realNormal[];

out vec3 ginitialVertPos;
out vec3 grealNormal;
out vec4 gcolor;


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

void sendInfoVertex(vec4 vecPos) {
    float displacementStrength = 1.0;
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    grealNormal = getNormal(vecPos.xy);
    gcolor = vec4(texture(heightmapFieldTex, vecPos.xy/texSize).rgb, 1.0);
    vecPos += vec4(grealNormal * getDisplacementLength(vecPos.xy) * displacementStrength, 0.0);
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
    float displacementStrength = 1.0;
    vec2 texSize = textureSize(heightmapFieldTex, 0);

    vec2 vecPos = initialVertPos[0].xy;
    if (vecPos.x >= texSize.x - 1 || vecPos.y >= texSize.y - 1) return;
    if (vecPos.x < min_vertice_positions.x || vecPos.x > max_vertice_positions.x ||
        vecPos.y < min_vertice_positions.y || vecPos.y > max_vertice_positions.y) return;

//    subdivision(vecPos);
//    return;

    vec4 v1 = vec4(vecPos + vec2(0, 0), getHeight(vecPos + vec2(0, 0)), 1.0);
    vec4 v2 = vec4(vecPos + vec2(1, 0), getHeight(vecPos + vec2(1, 0)), 1.0);
    vec4 v3 = vec4(vecPos + vec2(0, 1), getHeight(vecPos + vec2(0, 1)), 1.0);
    vec4 v4 = vec4(vecPos + vec2(1, 1), getHeight(vecPos + vec2(1, 1)), 1.0);

//    grealNormal = normalize(cross(v1.xyz - v2.xyz, v3.xyz - v1.xyz)); //realNormal[0];

    grealNormal = getNormal(v1.xy);
    v1 += vec4(grealNormal * getDisplacementLength(v1.xy) * displacementStrength, 0.0);
    ginitialVertPos = v1.xyz;
    gcolor = vec4(texture(heightmapFieldTex, v1.xy/texSize).rgb, 1.0);
    gl_Position = proj_matrix * mv_matrix * v1;
    EmitVertex();

    grealNormal = getNormal(v2.xy);
    v2 += vec4(grealNormal * getDisplacementLength(v2.xy) * displacementStrength, 0.0);
    ginitialVertPos = v2.xyz;
    gcolor = vec4(texture(heightmapFieldTex, v2.xy/texSize).rgb, 1.0);
    gl_Position = proj_matrix * mv_matrix * v2;
    EmitVertex();

    grealNormal = getNormal(v3.xy);
    v3 += vec4(grealNormal * getDisplacementLength(v3.xy) * displacementStrength, 0.0);
    ginitialVertPos = v3.xyz;
    gcolor = vec4(texture(heightmapFieldTex, v3.xy/texSize).rgb, 1.0);
    gl_Position = proj_matrix * mv_matrix * v3;
    EmitVertex();

    grealNormal = getNormal(v4.xy);
    v4 += vec4(grealNormal * getDisplacementLength(v4.xy) * displacementStrength, 0.0);
    ginitialVertPos = v4.xyz;
    gcolor = vec4(texture(heightmapFieldTex, v4.xy/texSize).rgb, 1.0);
    gl_Position = proj_matrix * mv_matrix * v4;
    EmitVertex();
    EndPrimitive();
}
