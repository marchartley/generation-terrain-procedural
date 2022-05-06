#version 120
// layout (location=0) in vec3 position;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

void main(void)
{
    vec3 position = gl_Vertex.xyz;
    gl_Position = proj_matrix * mv_matrix * vec4(position, 1.0);
}
