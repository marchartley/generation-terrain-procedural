#version 430
layout (location=0) in vec3 position;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

out vec3 initialVertPos;

void main(void)
{
    initialVertPos = position;
    gl_Position = proj_matrix * mv_matrix * vec4(initialVertPos, 1.0);
}
