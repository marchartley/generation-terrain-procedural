#version 430

out vec4 fragColor;

in vec3 particleCenterPos;
uniform vec3 min_vertice_positions;
uniform vec3 max_vertice_positions;
uniform bool clipPlaneActive;
uniform vec3 clipPlanePosition;
uniform vec3 clipPlaneDirection;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec4 color = vec4(1.0, 0.0, 0.0, 1.0);

void main(void)
{

    if (min_vertice_positions.x > particleCenterPos.x || particleCenterPos.x > max_vertice_positions.x || min_vertice_positions.y > particleCenterPos.y || particleCenterPos.y > max_vertice_positions.y || min_vertice_positions.z > particleCenterPos.z || particleCenterPos.z > max_vertice_positions.z)
        discard;

    if (clipPlaneActive) {
        if (dot((particleCenterPos.xyz - clipPlanePosition), clipPlaneDirection) > 0) {
            discard;
        }
    }

    fragColor = color;

}
