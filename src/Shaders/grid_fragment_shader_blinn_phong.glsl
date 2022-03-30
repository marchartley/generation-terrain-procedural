#version 430
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

in vec4 varyingColor;
in vec3 varyingHalfH;

in vec3 varyingNormal;
in vec3 varyingLightDir;
in vec3 varyingVertPos;
in vec3 vertEyeSpacePos;
out vec4 fragColor;

in vec3 realNormal;

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
uniform Material ground_material;
uniform Material grass_material;
uniform float fogNear;
uniform float fogFar;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

void main(void)
{
    Material material = ground_material;
    vec3 N = normalize(varyingNormal);
    vec3 L = normalize(varyingLightDir);
    vec3 V = normalize(-varyingVertPos);
    vec3 R = reflect(N, -L);
    vec3 H = normalize(varyingHalfH);
    float cosTheta = dot(L, N);
    float cosPhi = dot(H, N);

//    float upValue = 1 - clamp(dot(3 * normalize(realNormal + .1 * vec3(random(varyingVertPos.yz), random(varyingVertPos.xz), random(varyingVertPos.xy))), normalize(vec3(0.0, 1.0, 0.0))), 0.0, 1.0);
    float upValue = pow(clamp(dot(normalize(realNormal), normalize(vec3(0.0, 0.0, 1.0))), 0.0, 1.0), 2);
    vec4 material_ambiant = (ground_material.ambiant * (1 - upValue) + grass_material.ambiant * upValue) * .8;
    vec4 material_diffuse = (ground_material.diffuse * (1 - upValue) + grass_material.diffuse * upValue) * .8;
    vec4 material_specular = (ground_material.specular * (1 - upValue) + grass_material.specular * upValue) * 1.;
    float material_shininness = ground_material.shininness * (1 - upValue) + grass_material.shininness * upValue;

    vec4 ambiant = ((globalAmbiant * material_ambiant) + (light.ambiant * material_ambiant));
    vec4 diffuse = light.diffuse * material_diffuse * max(cosTheta, 0.0);
//    vec4 ambiantColor = vec4(varyingColor.xyz * .5, 1.0);
//    vec4 ambiant = ((globalAmbiant * ambiantColor) + (light.ambiant * ambiantColor));
//    vec4 diffuse = light.diffuse * (varyingColor*1.1) * max(cosTheta, 0.0);
    vec4 specular = light.specular * material_specular * pow(max(cosPhi, 0.0), material_shininness * 3.0);

    vec4 material_color = vec4((ambiant + diffuse + specular).xyz, 1.0);

    vec4 fogColor = vec4(0.7, 0.8, 0.9, 1.0);
    float dist = length(vertEyeSpacePos);
    float fogFactor = clamp(((fogFar - dist) / (fogFar - fogNear)), 0.0, 1.0);

    fragColor = mix(fogColor, material_color, fogFactor);
//    fragColor = vec4(1.0, 1.0, 1.0, 1.0);

//    bool isLight = length(vec3(proj_matrix * vec4(light.position, 1.0)).yz*10 - vertEyeSpacePos.xy) < 1.0;
//    if (isLight)
//        fragColor = vec4(V, 1.0) ;//vec4(1.0, 1.0, 1.0, 1.0);

}
