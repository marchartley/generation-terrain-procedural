#version 430

varying vec4 fragColor;

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
uniform PositionalLight lights[100];
uniform int lights_count;
uniform Material material;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

in vec4 varyingColor;
in vec3 varyingHalfH;

in vec3 varyingNormal;
in vec3 varyingLightDir;
in vec3 varyingVertPos;
in vec3 vertEyeSpacePos;
in vec3 initialVertPos;

void main(void)
{
    fragColor = vec4(1.0, 1.0, 1.0, 1.0) * varyingColor;

    vec3 N = normalize(varyingNormal);
    vec3 L = normalize(varyingLightDir);
    vec3 V = normalize(-varyingVertPos);
    vec3 R = reflect(-L, N);
    vec3 H = normalize(varyingHalfH);
    float cosTheta = dot(L, N);
    float cosPhi = dot(H, N);

    vec4 ambiant = ((globalAmbiant * material.ambiant) + (light.ambiant * material.ambiant));

    vec4 diffuse = light.diffuse * material.diffuse * max(cosTheta, 0.0);
    vec4 specular = light.specular * material.specular * pow(max(cosPhi, 0.0), material.shininness);

    for (int iLight = 0; iLight < lights_count; iLight++) {
        vec3 iL = normalize(vec4(mv_matrix * vec4(lights[iLight].position.xyz - initialVertPos, 1.0)).xyz);
        vec3 iR = reflect(-iL, N);
        float iCosTheta = dot(iL, N);
        float iCosPhi = dot(iR, V);
        ambiant +=  (lights[iLight].ambiant * material.ambiant);
        diffuse += lights[iLight].diffuse * material.diffuse * max(iCosTheta, 0.0);
        specular += lights[iLight].specular * material.specular * pow(max(iCosPhi, 0.0), material.shininness);
    }

    fragColor = vec4(ambiant + diffuse + specular);

}
