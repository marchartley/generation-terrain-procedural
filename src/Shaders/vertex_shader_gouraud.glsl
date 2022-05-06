#version 430
layout (location=0) in vec3 position;
//layout (location=1) in vec3 texture;
layout (location=2) in vec3 normal;

struct PositionalLight {
    vec4 ambiant;
    vec4 diffuse;
    vec4 specular;
    vec3 position;
};

struct Material {
    vec4 ambiant;
    vec4 diffuse;
    vec4 specular;
    float shininness;
};

uniform vec4 globalAmbiant;
uniform PositionalLight light;
uniform Material material;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

out vec4 varyingColor;
out vec3 initialVertPos;
out vec3 varyingVertPos;
out vec3 varyingLightDir;
out vec3 varyingNormal;
out vec3 varyingHalfH;


void main(void)
{
    vec4 P = mv_matrix * vec4(position, 1.0);
    vec3 N = normalize(vec4(norm_matrix * vec4(normal, 1.0)).xyz);
    vec3 L = normalize(light.position - P.xyz);
    vec3 V = normalize(-P.xyz);
    vec3 R = reflect(-L, N);

    varyingColor = vec4(1.0, 1.0, 1.0, 1.0);

    vec3 light_position = vec4(mv_matrix * vec4(light.position, 1.0)).xyz;
    initialVertPos = vec3(position);
    varyingVertPos = vec4(mv_matrix * vec4(position, 1.0)).xyz;
    varyingLightDir = light_position - varyingVertPos;
    varyingNormal = vec4(transpose(inverse(mv_matrix)) * vec4(normal, 1.0)).xyz;
    varyingHalfH = (varyingLightDir - varyingVertPos).xyz;

    gl_Position = proj_matrix * mv_matrix * vec4(position, 1.0);

}
