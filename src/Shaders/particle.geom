#version 430

#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable

layout ( points ) in;
layout ( triangle_strip, max_vertices = 4 ) out;

in vec3 initialVertPos[];

//Geometry Shader entry point
void main(void) {
    float scale = 0.05; //30.0 * 1.0/length(initialVertPos[0]);
    gl_Position = gl_in[0].gl_Position + vec4(-scale, -scale, 0.0, 0.0);
    EmitVertex();
    gl_Position = gl_in[0].gl_Position + vec4(-scale,  scale, 0.0, 0.0);
    EmitVertex();
    gl_Position = gl_in[0].gl_Position + vec4( scale, -scale, 0.0, 0.0);
    EmitVertex();
    gl_Position = gl_in[0].gl_Position + vec4( scale,  scale, 0.0, 0.0);
    EmitVertex();
    EndPrimitive();
}
