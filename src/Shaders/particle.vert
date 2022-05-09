#version 430
layout (location=0) in vec3 position;
layout (location=3) in vec3 offset;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

void main(void)
{
    vec3 offsetPos = position + offset;
    gl_Position = proj_matrix * mv_matrix * vec4(offsetPos, 1.0);
}
