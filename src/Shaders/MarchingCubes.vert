#version 430
layout (location=0) in vec3 position;
//layout (location=1) in vec3 texture;
layout (location=2) in vec3 normal;
layout (location=3) in vec3 color;
layout (location=4) in vec3 instanceOffset;

float random (in vec2 st) {
    return fract(sin(dot(st.xy,
                         vec2(12.9898,78.233)))*
        43758.5453123);
}

// Based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float noise (in vec2 st) {
    vec2 i = floor(st);
    vec2 f = fract(st);

    // Four corners in 2D of a tile
    float a = random(i);
    float b = random(i + vec2(1.0, 0.0));
    float c = random(i + vec2(0.0, 1.0));
    float d = random(i + vec2(1.0, 1.0));

    vec2 u = f * f * (3.0 - 2.0 * f);

    return mix(a, b, u.x) +
            (c - a)* u.y * (1.0 - u.x) +
            (d - b) * u.x * u.y;
}

float fbm (in vec2 st) {
    // Initial values
    float value = 0.0;
    float amplitud = .5;
    float frequency = 0.;
    //
    // Loop of octaves
    for (int i = 0; i < 6; i++) {
        value += amplitud * noise(st);
        st *= 2.;
        amplitud *= .5;
    }
    return value;
}
float fbm3 (in vec3 st) {
    return fbm(st.xy + st.z);
}

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

uniform float fogNear;
uniform float fogFar;
uniform float offsetX;
uniform float offsetY;

uniform vec4 globalAmbiant;
uniform PositionalLight light;
uniform bool isSpotlight;
uniform Material ground_material;
uniform Material grass_material;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;


out vec3 initialVertPos;
out vec3 realNormal;

void main(void)
{
    vec3 position = vec3(position.x + offsetX + instanceOffset.x, position.y + offsetY + instanceOffset.y, position.z + instanceOffset.z);

    initialVertPos = vec3(position);
    realNormal = normal;

    gl_Position = proj_matrix * mv_matrix * vec4(position, 1.0);

}
