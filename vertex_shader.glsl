#version 430
layout (location=0) in vec3 position;
//layout (location=1) in vec3 texture;
//layout (location=3) in vec3 normal;

//struct PositionalLight {
//    vec4 ambiant;
//    vec4 diffuse;
//    vec4 specular;
//    vec3 position;
//};

//struct Material {
//    vec4 ambiant;
//    vec4 diffuse;
//    vec4 specular;
//    float shininness;
//};

uniform float offsetX;
uniform float offsetY;
//in float waterIndex = 30.0;

//uniform vec4 globalAmbiant;
//uniform PositionalLight light;
//uniform Material material;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

out vec4 varyingColor;

void main(void)
{
    float waterIndex = 20.0;
    vec3 position = vec3(position.x + offsetX, position.y + offsetY, position.z);
//    position = vec3(position.x, position.y, max(position.z, waterIndex-0.01));
//    vec4 color;
//    vec4 P = mv_matrix * vec4(position, 1.0);
//    vec3 N = normalize((norm_matrix * vec4(normal, 1.0)).xyz);

    gl_Position = proj_matrix * mv_matrix * vec4(position, 1.0);
    if (position.z < waterIndex) {
        varyingColor = vec4(0.0, 0.0, 1.0, 1.0);
//        gl_Position.z = waterIndex;
    } else {
        varyingColor = vec4(0.0, 1.0, 0.0, 1.0);
    }

}
