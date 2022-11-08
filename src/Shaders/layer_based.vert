#version 430
layout (location=0) in vec3 position;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform sampler3D matIndicesTex;
uniform sampler3D matHeightsTex;

out vec3 initialVertPos;

void main(void)
{
    initialVertPos = vec3(position);

    // Start from the ground (z = 0)
    vec3 ground_position = vec3(position.x, position.y, 0.0);


    gl_Position = proj_matrix * mv_matrix * vec4(ground_position, 1.0);
}
