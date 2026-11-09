#version 430

out vec4 fragColor;

in vec3 particleCenterPos;
in vec3 particleColor;
in vec2 uv;
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

//    if (min_vertice_positions.x > particleCenterPos.x || particleCenterPos.x > max_vertice_positions.x || min_vertice_positions.y > particleCenterPos.y || particleCenterPos.y > max_vertice_positions.y || min_vertice_positions.z > particleCenterPos.z || particleCenterPos.z > max_vertice_positions.z)
//        discard;

//    if (clipPlaneActive) {
//        if (dot((particleCenterPos.xyz - clipPlanePosition), clipPlaneDirection) > 0) {
//            discard;
//        }
//    }

    vec2 radial = uv * 2.0 - 1.0;
    float radiusSquared = dot(radial, radial);
    if (radiusSquared >= 1.0) discard;
    float alpha = 1.0 - sqrt(radiusSquared);
    fragColor = vec4(particleColor.xyz, alpha);

}
