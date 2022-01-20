#version 430
layout (location=0) in vec3 position;
//layout (location=1) in vec3 texture;
layout (location=2) in vec3 normal;
layout (location=3) in vec3 color;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

void main(void)
{
    gl_Position = proj_matrix * mv_matrix * vec4(position, 1.0);

}
