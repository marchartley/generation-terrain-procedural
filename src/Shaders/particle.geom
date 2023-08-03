#version 430

//#extension GL_EXT_geometry_shader4 : enable
//#extension GL_EXT_gpu_shader4 : enable

layout ( points ) in;
layout ( triangle_strip, max_vertices = 4 ) out;

in vec3 initialVertPos[];
in vec3 initialVertCol[];
out vec3 particleCenterPos;
out vec3 particleColor;
out vec2 uv;

//Geometry Shader entry point
void main(void) {
    float scale = 2.0; // 0.05; //30.0 * 1.0/length(initialVertPos[0]);
    particleCenterPos = initialVertPos[0].xyz;
    particleColor = initialVertCol[0];
    uv = vec2(0, 0);
    gl_Position = gl_in[0].gl_Position + vec4(-scale, -scale, 0.0, 0.0);
    EmitVertex();
    uv = vec2(0, 1);
    gl_Position = gl_in[0].gl_Position + vec4(-scale,  scale, 0.0, 0.0);
    EmitVertex();
    uv = vec2(1, 0);
    gl_Position = gl_in[0].gl_Position + vec4( scale, -scale, 0.0, 0.0);
    EmitVertex();
    uv = vec2(1, 1);
    gl_Position = gl_in[0].gl_Position + vec4( scale,  scale, 0.0, 0.0);
    EmitVertex();
    EndPrimitive();
}
