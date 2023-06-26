#version 430
layout (location=0) in vec3 position;
layout (location=3) in vec3 color;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform float time;
uniform float maxTerrainHeight;

out vec3 initialVertPos;
out vec3 initialVertCol;

void main(void)
{
    initialVertPos = position; // + vec3(0, 0, -0.05) * (time);
//    initialVertPos.z = mod(initialVertPos.z, maxTerrainHeight);
    initialVertCol = color;
    gl_Position = proj_matrix * mv_matrix * vec4(initialVertPos, 1.0);
}
