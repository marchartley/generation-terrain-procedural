#version 430
layout (location=0) in vec3 position;
layout (location=3) in vec3 colors;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

out vec3 vColor;

void main(void)
{
    gl_Position = proj_matrix * mv_matrix * vec4(position, 1.0);
    vColor = colors;
}
