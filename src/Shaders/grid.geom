#version 430
//New G80 extensions
#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable

layout ( points ) in;
layout ( triangle_strip, max_vertices = 4) out;

//height data field texture
uniform sampler2D heightmapFieldTex;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

//Vertices position for fragment shader
in vec3 initialVertPos[];
in vec3 realNormal[];

out vec3 ginitialVertPos;
out vec3 grealNormal;


float getHeight(vec2 pos) {
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    return (pos.x >= texSize.x || pos.y >= texSize.y) ? 0 : texture(heightmapFieldTex, pos/texSize).r;
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

//Geometry Shader entry point
void main(void) {
    vec2 texSize = textureSize(heightmapFieldTex, 0);

    vec2 vecPos = initialVertPos[0].xy;
    if (vecPos.x >= texSize.x - 1 || vecPos.y >= texSize.y - 1) return;

    vec4 v1 = vec4(vecPos + vec2(0, 0), getHeight(vecPos + vec2(0, 0)), 1.0);
    vec4 v2 = vec4(vecPos + vec2(1, 0), getHeight(vecPos + vec2(1, 0)), 1.0);
    vec4 v3 = vec4(vecPos + vec2(0, 1), getHeight(vecPos + vec2(0, 1)), 1.0);
    vec4 v4 = vec4(vecPos + vec2(1, 1), getHeight(vecPos + vec2(1, 1)), 1.0);

//    grealNormal = normalize(cross(v1.xyz - v2.xyz, v3.xyz - v1.xyz)); //realNormal[0];

    gl_Position = proj_matrix * mv_matrix * v1;
    grealNormal = getNormal(v1.xy);
    ginitialVertPos = v1.xyz;
    EmitVertex();

    gl_Position = proj_matrix * mv_matrix * v2;
    grealNormal = getNormal(v2.xy);
    ginitialVertPos = v2.xyz;
    EmitVertex();

    gl_Position = proj_matrix * mv_matrix * v3;
    grealNormal = getNormal(v3.xy);
    ginitialVertPos = v3.xyz;
    EmitVertex();

    gl_Position = proj_matrix * mv_matrix * v4;
    grealNormal = getNormal(v4.xy);
    ginitialVertPos = v4.xyz;
    EmitVertex();
    EndPrimitive();
}
