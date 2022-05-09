#version 430

out vec4 fragColor;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec3 offset;

uniform vec4 color = vec4(1.0, 1.0, 1.0, 1.0);

void main(void)
{
    fragColor = color;

}
