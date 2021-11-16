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

uniform float offsetX;
uniform float offsetY;
//in float waterIndex = 30.0;

uniform vec4 globalAmbiant;
uniform PositionalLight light;
uniform Material material;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

out vec4 varyingColor;


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



void main(void)
{
    float material_shininness = 0.01;
    vec4 groundColor = vec4(0.163, 0.106, 0.034, 1.0);
    vec4 grassColor = vec4(0.018, 0.447, 0.038, 1.0);

    vec4 fakePos = proj_matrix * mv_matrix * vec4(position, 1.0);

    vec3 upVector = vec3(0.0, 0.0, 1.0);
    float upVal = dot(normalize(normal), upVector);

    vec4 currentColor = mix(groundColor, grassColor, -.5 + upVal);
//    material.ambiant = vec4(mix(groundColor, grassColor, upVal), 1.0);
//    float waterIndex = 20.0;
    vec3 position = vec3(position.x + offsetX, position.y + offsetY, position.z);
//    position = vec3(position.x, position.y, max(position.z, waterIndex-0.01));
    vec4 color;
    vec4 P = mv_matrix * vec4(position, 1.0);
    vec3 N = normalize(vec4(norm_matrix * vec4(normal, 1.0)).xyz);
//    vec3 N = normal;
    vec3 L = normalize(light.position - P.xyz);
    vec3 V = normalize(-P.xyz);
    vec3 R = reflect(-L, N);
//    L = normalize(L + normalize(vec3(fbm(fakePos.yz / 10.0), fbm(fakePos.xz / 10.0), fbm(fakePos.xy / 10.0))));

//    vec3 ambiant = ((globalAmbiant * material.ambiant) + (light.ambiant * material.ambiant)).xyz;
    vec3 ambiant = (light.ambiant * currentColor).xyz;
//    vec3 diffuse = light.diffuse.xyz * material.diffuse.xyz * max(dot(N, L), 0.0);
    vec3 diffuse = (light.diffuse * currentColor).xyz * max(dot(N, L), 0.0);
    vec3 specular = light.specular.xyz * material.specular.xyz;
//    vec3 specular = vec3(0.0);
    if (dot(R, V) > 0.0) {
        specular = specular * pow(dot(R, V), material.shininness / 10.0); //material.shininness);
    }
    else {
        specular = specular * 0.0;
    }
//    * pow(max(dot(R, V), 0.0), material.shininness);

    ambiant *= fbm3(position.xyz);
//    diffuse *= fbm3(position.xyz);
//    specular *= fbm3(position.xyz);

    vec3 tmp_color = vec3(0.0);
    tmp_color = ambiant;
    tmp_color = tmp_color + diffuse;
    tmp_color = tmp_color + specular;
//    tmp_color = tmp_color + abs(specular) / 100.0;
    varyingColor = vec4(tmp_color, 1.0);
//    varyingColor = vec4((ambiant + diffuse + specular), 1.0);
//    varyingColor = vec4(specular, 1.0);

    gl_Position = proj_matrix * mv_matrix * vec4(position, 1.0);
//    if (position.z < waterIndex) {
//        varyingColor = vec4(0.0, 0.0, 1.0, 1.0);
////        gl_Position.z = waterIndex;
//    } else {
//        varyingColor = vec4(0.0, 1.0, 0.0, 1.0);
//    }

}
